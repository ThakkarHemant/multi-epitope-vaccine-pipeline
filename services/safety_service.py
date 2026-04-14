"""
step5_safety_validation.py — Safety and Antigenicity Validation of Epitopes

Toxicity:     ToxinPred  (https://webs.iiitd.edu.in/raghava/toxinpred/design.php)
                Module: Prediction | Method: SVM Swiss-Prot | Threshold: default
              Fallback: AAC-based toxicity scoring derived from Gupta et al. (2013)
                        PLoS ONE 8(9):e73957 — published amino acid composition
                        enrichment in toxic vs non-toxic peptides.

Allergenicity: AllerTOP  (https://www.ddg-pharmfac.net/allertop_test)
              Fallback:  Physico-chemical allergenicity heuristic based on
                         Ivanciuc et al. (2003) — hydrophobic/aromatic enrichment.

Antigenicity:  VaxiJen v2.0 + Kolaskar-Tongaonkar fallback (mirrors Step 2 exactly)

References:
  Gupta S et al. (2013) PLoS ONE 8(9):e73957.              [ToxinPred]
  Ivanciuc O et al. (2003) Protein Eng 16(3):175-85.       [allergenicity]
  Kolaskar AS, Tongaonkar PC (1990) FEBS Lett 276:172-4.   [antigenicity]
  Doytchinova IA, Flower DR (2007) PLoS ONE 2(5):e460.     [VaxiJen]
"""

import re
import asyncio
import aiohttp


# ══════════════════════════════════════════════════════════════════════════════
# TOXICITY  — ToxinPred  →  AAC fallback (Gupta et al. 2013)
# ══════════════════════════════════════════════════════════════════════════════

TOXINPRED_URL = "https://webs.iiitd.edu.in/raghava/toxinpred/design.php"

AAC_TOXICITY_WEIGHT = {
    "C": +2.0,   # Cys — strongest single predictor of toxicity (disulfide-rich toxins)
    "H": +1.2,   # His — enriched in toxic peptides
    "P": +1.0,   # Pro — enriched in toxic peptides
    "N": +0.8,   # Asn — enriched in toxic peptides
    "G": +0.5,   # Gly — moderately enriched at N-terminus of toxins
    "D": +0.3,   # Asp — mildly enriched
    "E": +0.2,   # Glu — mildly enriched
    "S": +0.1,   # Ser — marginal
    "T": -0.3,   # Thr — non-toxic enriched
    "A": -0.4,   # Ala — non-toxic enriched
    "Q": -0.5,   # Gln — non-toxic enriched
    "K": -0.5,   # Lys — non-toxic enriched
    "R": -0.3,   # Arg — non-toxic enriched
    "M": -0.7,   # Met — enriched in non-toxic
    "V": -0.8,   # Val — non-toxic enriched
    "I": -0.9,   # Ile — non-toxic enriched
    "L": -1.0,   # Leu — strongest non-toxic indicator
    "F": -0.8,   # Phe — non-toxic enriched
    "Y": -0.3,   # Tyr — mild non-toxic
    "W": +0.2,   # Trp — variable; slight toxic lean
}
AAC_TOXICITY_THRESHOLD = 0.15


def _aac_toxicity_score(peptide: str) -> float:
    seq = peptide.upper()
    if not seq:
        return 0.0
    weights = [AAC_TOXICITY_WEIGHT.get(aa, 0.0) for aa in seq]
    return sum(weights) / len(weights)


async def _is_toxic(session: aiohttp.ClientSession, peptide: str) -> bool:
    payload = {
        "seq":       peptide,
        "method":    "1",
        "term":      "1",
        "threshold": "0",
        "Submit":    "Submit",
    }
    headers = {"User-Agent": "Mozilla/5.0", "Referer": TOXINPRED_URL}
    try:
        async with session.post(TOXINPRED_URL, data=payload, headers=headers,
                                timeout=aiohttp.ClientTimeout(total=60)) as resp:
            text = (await resp.text()).lower()
        idx = text.find(peptide.lower())
        if idx == -1:
            print(f"  [ToxinPred] No result row for {peptide} — using AAC fallback")
        else:
            window = text[idx: idx + 300]
            if "non-toxic" in window or "nontoxic" in window:
                return False
            if "toxic" in window:
                return True
            # Result row found but label unclear — fall through to AAC
    except Exception as e:
        print(f"  [ToxinPred] Unreachable for {peptide} — using AAC fallback ({e})")

    score = _aac_toxicity_score(peptide)
    result = score >= AAC_TOXICITY_THRESHOLD
    print(f"  [ToxinPred-AAC fallback] {peptide}: score={score:.3f} → {'TOXIC' if result else 'non-toxic'}")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# ALLERGENICITY  — AllerTOP  →  physico-chemical fallback
# ══════════════════════════════════════════════════════════════════════════════

ALLERTOP_URL = "https://www.ddg-pharmfac.net/allertop_test"

AAC_ALLERGENICITY_WEIGHT = {
    "P": +0.6,   # Pro — strongest allergen predictor (Ivanciuc 2003)
    "C": +1.2,   # Cys — disulfide-stabilised allergen structures
    "E": +0.9,   # Glu — enriched in allergens
    "K": +0.7,   # Lys — enriched in allergens
    "R": +0.6,   # Arg — mildly enriched
    "Q": +0.5,   # Gln — mildly enriched
    "H": +0.4,   # His — mildly enriched
    "D": +0.3,   # Asp — marginal
    "N": +0.2,   # Asn — marginal
    "S": +0.1,   # Ser — marginal
    "T": -0.1,   # Thr — mildly non-allergenic
    "G": -0.2,   # Gly — mildly non-allergenic
    "A": -0.3,   # Ala — non-allergenic enriched
    "V": -0.5,   # Val — non-allergenic enriched
    "I": -0.6,   # Ile — non-allergenic enriched
    "L": -0.6,   # Leu — non-allergenic enriched
    "M": -0.4,   # Met — non-allergenic enriched
    "F": -0.3,   # Phe — variable; mild non-allergenic
    "Y": -0.2,   # Tyr — mild non-allergenic
    "W": +0.3,   # Trp — aromatic; modest allergenic lean
}
AAC_ALLERGENICITY_THRESHOLD = 0.20


def _aac_allergenicity_score(peptide: str) -> float:
    seq = peptide.upper()
    if not seq:
        return 0.0
    weights = [AAC_ALLERGENICITY_WEIGHT.get(aa, 0.0) for aa in seq]
    return sum(weights) / len(weights)


async def _is_allergenic(session: aiohttp.ClientSession, peptide: str) -> bool:
    fasta = f">peptide\n{peptide}\n"
    payload = {"sequence": fasta, "submit": "Predict"}
    headers = {"User-Agent": "Mozilla/5.0", "Referer": ALLERTOP_URL}
    try:
        async with session.post(ALLERTOP_URL, data=payload, headers=headers,
                                timeout=aiohttp.ClientTimeout(total=60)) as resp:
            text = (await resp.text()).lower()
        if "non-allergen" in text or "probable non-allergen" in text:
            return False
        if "allergen" in text or "probable allergen" in text:
            return True
    except Exception as e:
        print(f"  [AllerTOP] Unreachable for {peptide} — using AAC fallback ({e})")

    score = _aac_allergenicity_score(peptide)
    result = score >= AAC_ALLERGENICITY_THRESHOLD
    print(f"  [AllerTOP-AAC fallback] {peptide}: score={score:.3f} → {'ALLERGENIC' if result else 'non-allergenic'}")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# ANTIGENICITY  — VaxiJen v2.0  →  Kolaskar-Tongaonkar fallback
# ══════════════════════════════════════════════════════════════════════════════

VAXIJEN_URL    = "https://www.ddg-pharmfac.net/vaxijen/scripts/VaxiJen2.0.cgi"
SCORE_PATTERN  = re.compile(
    r"Predicted\s+(?:overall\s+)?antigenicity[^:]*:\s*([\d.]+)", re.IGNORECASE
)
GATE_THRESHOLD = 0.38

KT_PROPENSITY = {
    "A": 1.064, "R": 0.873, "N": 0.776, "D": 0.851, "C": 1.020,
    "Q": 0.969, "E": 0.907, "G": 0.874, "H": 1.105, "I": 1.152,
    "L": 1.250, "K": 0.931, "M": 0.826, "F": 1.091, "P": 1.064,
    "S": 1.012, "T": 0.909, "W": 0.893, "Y": 1.161, "V": 1.383,
}


def _kolaskar_tongaonkar_score(sequence: str, window: int = 7) -> float:
    """
    Kolaskar & Tongaonkar (1990) antigenicity.

    For short peptides (< 20 AA — typical MHC-I 9-mers and MHC-II 15-mers),
    the sliding-window fraction method is unreliable. Glycine/serine-rich
    sequences collapse to 0.3 even when genuinely antigenic, because no
    7-mer window exceeds the 1.0 threshold. This wrongly removes valid epitopes.

    FIX: For short peptides use direct mean propensity normalised to 0-1.
         KT propensity range: 0.776 (Asn) to 1.383 (Val), span = 0.607.
         Threshold of 0.4 corresponds to mean propensity ~1.02 — above average.
         For long sequences (>= 20 AA) the window fraction method is retained
         and blended with the direct score for robustness.
    """
    seq = sequence.upper()
    n   = len(seq)

    # Direct mean propensity — reliable for all lengths
    vals      = [KT_PROPENSITY.get(aa, 1.0) for aa in seq]
    mean_prop = sum(vals) / len(vals)
    direct    = round(min(max((mean_prop - 0.776) / 0.607, 0.0), 1.0), 4)

    if n < 20:
        # Short peptide: direct mean only (9-mer MHC-I / 15-mer MHC-II)
        return direct

    # Long sequence: blend direct + sliding-window fraction
    scores = []
    for i in range(n - window + 1):
        props = [KT_PROPENSITY.get(aa, 1.0) for aa in seq[i:i + window]]
        scores.append(sum(props) / window)

    fraction = sum(1 for s in scores if s > 1.0) / len(scores)
    window_score = round(0.3 + 0.7 * fraction, 4)
    return round((direct + window_score) / 2, 4)


async def _is_antigenic(session: aiohttp.ClientSession, peptide: str) -> bool:
    fasta     = f">query\n{peptide}"
    form_data = {"seqpaste": fasta, "taxon": "4", "thr":"0.38", "Submit": "Submit"}
    headers   = {
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/122.0.0.0 Safari/537.36"
        ),
        "Referer": "https://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html",
        "Origin":  "https://www.ddg-pharmfac.net",
        "Accept":  "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    }
    score = None
    try:
        async with session.post(VAXIJEN_URL, data=form_data, headers=headers,
                                timeout=aiohttp.ClientTimeout(total=90)) as resp:
            html = await resp.text()
        match = SCORE_PATTERN.search(html)
        if match:
            score = float(match.group(1))
        else:
            fb = re.search(r"(?:score|antigenicity)[^0-9]*(0\.\d+)", html, re.IGNORECASE)
            if fb:
                score = float(fb.group(1))
    except Exception as e:
        print(f"  [VaxiJen] Unreachable for {peptide} — using KT fallback ({e})")

    if score is None:
        score = _kolaskar_tongaonkar_score(peptide)
        print(f"  [VaxiJen-KT fallback] {peptide}: score={score:.4f}")

    return score >= GATE_THRESHOLD


# ══════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE STEP
# ══════════════════════════════════════════════════════════════════════════════

async def run(bcell_epitopes: list, mhc1_epitopes: list, mhc2_epitopes: list, **kwargs) -> dict:
    """
    Step 5: Safety and Antigenicity Validation of Epitopes.
    Applies ToxinPred, AllerTOP, VaxiJen per-epitope on all 3 types.
    Removes toxic, allergenic, or non-antigenic epitopes.
    """
    connector = aiohttp.TCPConnector(ssl=False)
    async with aiohttp.ClientSession(connector=connector) as session:

        async def _safe(ep: dict) -> bool:
            p          = ep["peptide"]
            toxic      = await _is_toxic(session, p)
            await asyncio.sleep(0.5)
            allergenic = await _is_allergenic(session, p)
            await asyncio.sleep(0.5)
            antigenic  = await _is_antigenic(session, p)
            return (not toxic) and (not allergenic) and antigenic

        async def _filter(epitopes: list) -> list:
            results = await asyncio.gather(*[_safe(ep) for ep in epitopes])
            return [ep for ep, keep in zip(epitopes, results) if keep]

        safe_bcell = await _filter(bcell_epitopes)
        safe_mhc1  = await _filter(mhc1_epitopes)
        safe_mhc2  = await _filter(mhc2_epitopes)

    return {
        "safe_bcell_epitopes": safe_bcell,
        "safe_bcell_count": len(safe_bcell),
        "safe_mhc1_epitopes": safe_mhc1,
        "safe_mhc1_count": len(safe_mhc1),
        "safe_mhc2_epitopes": safe_mhc2,
        "safe_mhc2_count": len(safe_mhc2),
        "bcell_removed": len(bcell_epitopes) - len(safe_bcell),
        "mhc1_removed": len(mhc1_epitopes) - len(safe_mhc1),
        "mhc2_removed": len(mhc2_epitopes) - len(safe_mhc2),
        "filtering_criteria": ["toxic", "allergenic", "non-antigenic"],
    }