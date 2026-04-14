"""
ifn_gamma_service.py — Step 7: IFN-Gamma Epitope Prediction

Primary:   IFNepitope web server (IMTECH/IIITD) — Hybrid method (motif + SVM)
           https://webs.iiitd.edu.in/raghava/ifnepitope/run_submit-old.php

           The server is a JOB-QUEUE system:
             POST  → run_submit-old.php  returns a page containing a results URL
                     (typically /raghava/ifnepitope/tmp/<job_id>/result.html
                      or a meta-refresh / anchor pointing to it)
             GET   → poll that URL until the result table appears (up to ~120 s)

Fallback 1: ifnepitope2 standalone package (pip install ifnepitope2)
Fallback 2: AAC-based scoring (Dhanda et al. 2013 Biol Direct 8:30)

References:
  Dhanda SK et al. (2013) Biol Direct 8:30.
  Dhall A et al. (2024) Sci Rep 14(1):26859.
"""

import asyncio
import logging
import re

import httpx

logger = logging.getLogger(__name__)

# ── Server constants ──────────────────────────────────────────────────────────
IFNEPITOPE_BASE       = "https://webs.iiitd.edu.in"
IFNEPITOPE_SUBMIT_URL = f"{IFNEPITOPE_BASE}/raghava/ifnepitope/run_submit-old.php"

POLL_INTERVAL_S = 5
POLL_MAX_WAIT_S = 120
POLL_MAX_TRIES  = POLL_MAX_WAIT_S // POLL_INTERVAL_S

# Matches a result-table row: SEQUENCE  0.1234  IFN-gamma inducer
SCORE_ROW_RE = re.compile(
    r"([A-Z]{6,25})\s+([\-]?\d+\.\d+)\s+"
    r"(IFN[\-\s]gamma\s+\w+|[Nn]on[\-\s][Ii]nducer|\w*[Ii]nducer)",
    re.IGNORECASE,
)
JOB_RESULT_URL_RE = re.compile(
    r'(?:href|url|content)[=\s:"\']+([^"\'>\s]*(?:result|output)[^"\'>\s]*\.(?:html?|php))',
    re.IGNORECASE,
)
META_REFRESH_RE = re.compile(
    r'<meta[^>]+http-equiv=["\']refresh["\'][^>]+content=["\'][^;]+;\s*url=([^"\'>\s]+)',
    re.IGNORECASE,
)


# ══════════════════════════════════════════════════════════════════════════════
# FALLBACK 2 — AAC scoring (Dhanda et al. 2013)
# ══════════════════════════════════════════════════════════════════════════════

IFN_AAC_WEIGHT: dict[str, float] = {
    "K": +1.5, "L": +1.2, "I": +0.9, "R": +0.8, "N": +0.6,
    "P": +0.5, "F": +0.4, "V": +0.3, "M": +0.2, "W": +0.1,
    "A": -0.1, "T": -0.2, "S": -0.3, "G": -0.4, "Q": -0.5,
    "H": -0.3, "D": -0.6, "E": -0.8, "C": -0.3, "Y": -0.2,
}
IFN_AAC_THRESHOLD = 0.10


def _aac_ifn_score(peptide: str) -> float:
    seq = peptide.upper()
    if not seq:
        return 0.0
    return round(sum(IFN_AAC_WEIGHT.get(aa, 0.0) for aa in seq) / len(seq), 4)


# ══════════════════════════════════════════════════════════════════════════════
# FALLBACK 1 — ifnepitope2 local package
# ══════════════════════════════════════════════════════════════════════════════

def _try_ifnepitope2_local(peptides: list[str]) -> dict[str, tuple[float, str]] | None:
    try:
        from ifnepitope2 import IFNepitope2  # type: ignore
    except ImportError:
        return None
    try:
        model = IFNepitope2(host=1)
        results = {}
        for pep in peptides:
            score, label = model.predict(pep)
            results[pep] = (round(float(score), 4), label)
        logger.info("ifnepitope2 local: predicted %d peptides", len(results))
        return results
    except Exception as exc:
        logger.warning("ifnepitope2 local failed: %s", exc)
        return None


# ══════════════════════════════════════════════════════════════════════════════
# PRIMARY — IFNepitope web server (job-queue aware)
# ══════════════════════════════════════════════════════════════════════════════

def _build_fasta(peptide: str, index: int) -> str:
    return f">Epitope_{index}\n{peptide.upper()}"


def _extract_result_url(html: str) -> str | None:
    """Pull the result-page URL out of the job submission acknowledgement."""
    for pattern in (META_REFRESH_RE, JOB_RESULT_URL_RE):
        m = pattern.search(html)
        if m:
            url = m.group(1).strip()
            if url.startswith("http"):
                return url
            return IFNEPITOPE_BASE + ("" if url.startswith("/") else "/") + url

    # Fallback: any /tmp/ path in the response
    tmp = re.search(r'(/raghava/ifnepitope/tmp/[^"\'>\s]+)', html, re.IGNORECASE)
    if tmp:
        return IFNEPITOPE_BASE + tmp.group(1)

    return None


def _parse_result_html(html: str) -> tuple[float, str] | None:
    """
    Return (score, label) if the result table is present, else None.
    None → result page not ready yet; caller should keep polling.
    """
    m = SCORE_ROW_RE.search(html)
    if m:
        return float(m.group(2)), m.group(3).strip()

    m2 = re.search(
        r"(IFN[\-\s]gamma\s+inducer|[Nn]on[\-\s][Ii]nducer)"
        r"[^\d]*([\-]?\d+\.\d+)|([\-]?\d+\.\d+)[^\d]*"
        r"(IFN[\-\s]gamma\s+inducer|[Nn]on[\-\s][Ii]nducer)",
        html, re.IGNORECASE,
    )
    if m2:
        g = m2.groups()
        return float(g[1] or g[2] or 0), (g[0] or g[3] or "Unknown").strip()

    return None


async def _predict_single(
    peptide: str,
    index: int,
    client: httpx.AsyncClient,
) -> tuple[float, str]:
    """
    Submit one peptide, detect the job-queue result URL, poll until ready.
    """
    fasta_text = _build_fasta(peptide, index)

    # Step 1 — submit
    resp = await client.post(
        IFNEPITOPE_SUBMIT_URL,
        files={
            "sequence": (None, fasta_text,        "text/plain"),
            "method":   (None, "hybrid",           "text/plain"),
            "model":    (None, "main",             "text/plain"),
            "job_name": (None, f"Job_{index}",     "text/plain"),
        },
        timeout=60,
    )
    resp.raise_for_status()
    submit_html = resp.text

    # Step 2 — maybe results are already in the submit response
    immediate = _parse_result_html(submit_html)
    if immediate:
        logger.debug("IFNepitope: %s — immediate result on submit page", peptide)
        return immediate

    # Step 3 — find the result URL
    result_url = _extract_result_url(submit_html)
    if not result_url and str(resp.url) != IFNEPITOPE_SUBMIT_URL:
        result_url = str(resp.url)

    if not result_url:
        raise ValueError(
            f"Could not locate result URL in IFNepitope submission response "
            f"for '{peptide}'. Snippet: {submit_html[:400]!r}"
        )

    logger.debug("IFNepitope: %s — polling %s", peptide, result_url)

    # Step 4 — poll
    for attempt in range(1, POLL_MAX_TRIES + 1):
        await asyncio.sleep(POLL_INTERVAL_S)
        try:
            poll_resp = await client.get(result_url, timeout=30)
            poll_resp.raise_for_status()
            parsed = _parse_result_html(poll_resp.text)
            if parsed:
                logger.debug("IFNepitope: %s ready on poll attempt %d", peptide, attempt)
                return parsed
            logger.debug("IFNepitope: %s not ready yet (%d/%d)", peptide, attempt, POLL_MAX_TRIES)
        except httpx.HTTPStatusError as exc:
            if exc.response.status_code == 404:
                continue   # job not written to disk yet
            raise

    raise TimeoutError(
        f"IFNepitope result not available after {POLL_MAX_WAIT_S}s for '{peptide}'. "
        f"Last polled: {result_url}"
    )


# ══════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE STEP
# ══════════════════════════════════════════════════════════════════════════════

async def run(safe_mhc2_epitopes: list, **kwargs) -> dict:
    """
    Step 7: IFN-Gamma Epitope Prediction.

    Priority:
      1. IFNepitope web server  (job-queue, polling)
      2. ifnepitope2 local pkg  (pip install ifnepitope2)
      3. AAC composition score  (Dhanda et al. 2013 — last resort)

    Returns:
      approach, model,
      ifn_gamma_epitopes  — inducers only  (for downstream steps)
      all_predictions     — EVERY peptide with score + label + tool
      ifn_gamma_count,
      total_predicted,
      errors
    """
    errors: list[dict] = []

    valid_eps = [ep for ep in safe_mhc2_epitopes if 8 <= len(ep["peptide"]) <= 25]
    skipped   = len(safe_mhc2_epitopes) - len(valid_eps)
    if skipped:
        logger.warning("IFNepitope: skipped %d peptides outside 8–25 AA range", skipped)
    logger.info("IFNepitope: predicting for %d MHC-II epitopes …", len(valid_eps))

    # ── 1. Web server — one shared 120 s budget across ALL peptides ──────────
    #    If the budget expires (or any hard error occurs) every un-scored
    #    peptide is sent straight to the fallback chain; no per-peptide wait.
    WEB_TOTAL_TIMEOUT_S = 120

    web_results: dict[str, tuple[float, str]] = {}
    web_failed:  list[dict] = []

    async def _run_all_web():
        async with httpx.AsyncClient(timeout=60, follow_redirects=True) as client:
            for idx, ep in enumerate(valid_eps, start=1):
                pep = ep["peptide"]
                try:
                    score, label = await _predict_single(pep, idx, client)
                    web_results[pep] = (score, label)
                    logger.info("IFNepitope web ✓  %s → %.4f (%s)", pep, score, label)
                except Exception as exc:
                    logger.warning("IFNepitope web ✗  %s — %s", pep, exc)
                    web_failed.append(ep)
                    errors.append({"peptide": pep, "source": "web", "error": str(exc)})

    try:
        await asyncio.wait_for(_run_all_web(), timeout=WEB_TOTAL_TIMEOUT_S)
    except asyncio.TimeoutError:
        # Budget exhausted — mark every peptide not yet scored as failed
        scored = set(web_results)
        timed_out = [ep for ep in valid_eps if ep["peptide"] not in scored
                     and ep not in web_failed]
        for ep in timed_out:
            web_failed.append(ep)
            errors.append({
                "peptide": ep["peptide"],
                "source":  "web",
                "error":   f"Web server did not respond within {WEB_TOTAL_TIMEOUT_S}s — falling back",
            })
        logger.warning(
            "IFNepitope: 120 s budget exhausted — %d peptide(s) sent to fallback",
            len(timed_out),
        )

    # ── 2. ifnepitope2 local ──────────────────────────────────────────────────
    local_results: dict[str, tuple[float, str]] = {}
    still_failed:  list[dict] = []

    if web_failed:
        local_raw = _try_ifnepitope2_local([ep["peptide"] for ep in web_failed])
        if local_raw:
            local_results = local_raw
        else:
            still_failed = web_failed
            logger.warning(
                "ifnepitope2 not installed — using AAC fallback for %d peptides. "
                "Install with: pip install ifnepitope2",
                len(still_failed),
            )

    # ── 3. AAC fallback ───────────────────────────────────────────────────────
    aac_results: dict[str, tuple[float, str]] = {}
    for ep in still_failed:
        pep   = ep["peptide"]
        score = _aac_ifn_score(pep)
        label = "IFN-gamma inducer" if score >= IFN_AAC_THRESHOLD else "Non-inducer"
        aac_results[pep] = (score, label)
        logger.info("[AAC fallback] %s → score=%.4f  label=%s", pep, score, label)

    # ── Merge ─────────────────────────────────────────────────────────────────
    all_results = {**web_results, **local_results, **aac_results}
    tool_map: dict[str, str] = {
        **{p: "IFNepitope (Hybrid web)"     for p in web_results},
        **{p: "IFNepitope2 (local package)" for p in local_results},
        **{p: "AAC fallback (Dhanda 2013)"  for p in aac_results},
    }

    # ── Build output ──────────────────────────────────────────────────────────
    ifn_gamma_epitopes: list[dict] = []
    all_predictions:    list[dict] = []

    for ep in valid_eps:
        pep = ep["peptide"]
        if pep not in all_results:
            continue

        score, label = all_results[pep]
        is_inducer   = "inducer" in label.lower() and "non" not in label.lower()

        record = {
            **ep,
            "ifn_gamma_score": round(score, 4),
            "prediction":      label,
            "is_inducer":      is_inducer,
            "tool":            tool_map.get(pep, "Unknown"),
        }
        all_predictions.append(record)
        if is_inducer:
            ifn_gamma_epitopes.append(record)

    ifn_gamma_epitopes.sort(key=lambda x: -x.get("ifn_gamma_score", 0))
    all_predictions.sort(key=lambda x: -x.get("ifn_gamma_score", 0))

    logger.info(
        "IFNepitope: %d inducers / %d total predicted",
        len(ifn_gamma_epitopes), len(all_predictions),
    )

    return {
        "approach":           "Motif and SVM hybrid",
        "model":              "IFN-gamma versus Non-IFN-gamma",
        "ifn_gamma_epitopes": ifn_gamma_epitopes,   # inducers only — downstream
        "all_predictions":    all_predictions,       # every peptide — report
        "ifn_gamma_count":    len(ifn_gamma_epitopes),
        "total_predicted":    len(all_predictions),
        "errors":             errors,
    }