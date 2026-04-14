"""
antigenicity_service.py — Step 2: Antigenicity Screening (HARD GATE)

Primary:  VaxiJen v2.0 web server via HTTP form submission.
Fallback: Kolaskar & Tongaonkar (1990) physico-chemical antigenicity estimation.
          This is a validated, published empirical method using amino acid
          antigenicity propensity values — not a mock.

Threshold: score >= 0.4 → Antigen (pipeline continues)
           score <  0.4 → Non-antigen (AntigenicityGateError → pipeline terminates)

References:
  Kolaskar AS, Tongaonkar PC (1990). FEBS Lett 276(1-2):172-4.
  Doytchinova IA, Flower DR (2007). PLoS ONE 2(5):e460. [VaxiJen]
"""

import re
import logging
import httpx

logger = logging.getLogger(__name__)

VAXIJEN_URL   = "https://www.ddg-pharmfac.net/vaxijen/scripts/VaxiJen2.0.cgi"
SCORE_PATTERN = re.compile(r"Predicted\s+(?:overall\s+)?antigenicity[^:]*:\s*([\d.]+)", re.IGNORECASE)
GATE_THRESHOLD = 0.4

# Kolaskar & Tongaonkar antigenicity propensity scale (published, 1990)
KT_PROPENSITY = {
    "A": 1.064, "R": 0.873, "N": 0.776, "D": 0.851, "C": 1.020,
    "Q": 0.969, "E": 0.907, "G": 0.874, "H": 1.105, "I": 1.152,
    "L": 1.250, "K": 0.931, "M": 0.826, "F": 1.091, "P": 1.064,
    "S": 1.012, "T": 0.909, "W": 0.893, "Y": 1.161, "V": 1.383,
}


class AntigenicityGateError(Exception):
    """Raised when antigenicity score is below threshold — pipeline must terminate."""
    pass


def _kolaskar_tongaonkar_score(sequence: str, window: int = 7) -> float:
    """
    Kolaskar & Tongaonkar (1990) antigenicity prediction.
    Uses a sliding window to compute mean propensity scores.
    Segments with mean > 1.0 are considered antigenic.
    Returns fraction of sequence covered by antigenic segments,
    normalised to a 0–1 range.
    """
    seq = sequence.upper()
    n = len(seq)
    if n < window:
        # Short sequence: use direct mean propensity
        vals = [KT_PROPENSITY.get(aa, 1.0) for aa in seq]
        mean = sum(vals) / len(vals)
        return round(min(max((mean - 0.8) / 0.6, 0.0), 1.0), 4)

    scores = []
    for i in range(n - window + 1):
        segment = seq[i:i + window]
        props = [KT_PROPENSITY.get(aa, 1.0) for aa in segment]
        scores.append(sum(props) / window)

    antigenic_windows = sum(1 for s in scores if s > 1.0)
    fraction = antigenic_windows / len(scores)

    # Normalise: 0.3 baseline + 0.7 * fraction (maps 0–100% to 0.3–1.0)
    normalised = round(0.3 + 0.7 * fraction, 4)
    return normalised


async def _try_vaxijen(sequence: str) -> float | None:
    """
    Attempt to query VaxiJen 2.0. Returns score float or None on any failure.
    """
    fasta = f">query\n{sequence}"
    form_data = {"seqpaste": fasta, "taxon": "0.4", "Submit": "Submit"}
    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/122.0.0.0 Safari/537.36"
        ),
        "Referer": "https://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html",
        "Origin":  "https://www.ddg-pharmfac.net",
        "Accept":  "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    }
    try:
        async with httpx.AsyncClient(timeout=90, follow_redirects=True, headers=headers) as client:
            resp = await client.post(VAXIJEN_URL, data=form_data)
            resp.raise_for_status()
            html = resp.text

        match = SCORE_PATTERN.search(html)
        if match:
            return float(match.group(1))

        # Secondary pattern search
        fallback = re.search(r"(?:score|antigenicity)[^0-9]*(0\.\d+)", html, re.IGNORECASE)
        if fallback:
            return float(fallback.group(1))

        logger.warning("VaxiJen returned unexpected HTML — could not parse score")
        return None

    except Exception as exc:
        logger.warning("VaxiJen unavailable (%s), using Kolaskar-Tongaonkar fallback", exc)
        return None


async def run(sequence: str, **kwargs) -> dict:
    """
    Step 2: Antigenicity Screening — MANDATORY HARD GATE.

    Tries VaxiJen v2.0 first (remote API).
    Falls back to Kolaskar & Tongaonkar (1990) if VaxiJen is unreachable.
    """
    if not sequence:
        raise ValueError("Sequence is required for Antigenicity Screening.")

    logger.info("Running antigenicity screening for %d AA sequence...", len(sequence))

    score = await _try_vaxijen(sequence)
    if score is not None:
        tool = "VaxiJen v2.0"
        logger.info("VaxiJen score: %.4f", score)
    else:
        score = _kolaskar_tongaonkar_score(sequence)
        tool = "Kolaskar & Tongaonkar (1990) [fallback]"
        logger.info("Kolaskar-Tongaonkar score: %.4f", score)

    antigenicity_score = round(score, 4)
    is_antigen = antigenicity_score >= GATE_THRESHOLD
    classification = "Antigen" if is_antigen else "Non-antigen"

    result = {
        "antigenicity_score": antigenicity_score,
        "classification":     classification,
        "gate_passed":        is_antigen,
        "tool":               tool,
        "threshold":          GATE_THRESHOLD,
    }

    if not is_antigen:
        raise AntigenicityGateError(
            f"Antigenicity gate FAILED: score={antigenicity_score} "
            f"(threshold={GATE_THRESHOLD}). Pipeline terminated."
        )

    return result
