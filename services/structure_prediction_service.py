"""
structure_prediction_service.py — Structure Prediction

Primary:  ESMFold API (Meta AI)
          https://api.esmatlas.com/foldSequence/v1/pdb/
          - Free, no API key required
          - POST raw sequence → returns real PDB with pLDDT B-factors
          - Length limit: 400 AA (vaccine constructs are ~100-200 AA)
          - Confidence: mean pLDDT parsed from B-factor column (0-100 scale)

Fallback: Clearly-labelled placeholder PDB with confidence_score = None
          so the pipeline never crashes on timeout/error.

References:
  Lin Z et al. (2023) Science 379(6637):1123-30.  [ESMFold]
"""

import re
import httpx

ESMFOLD_URL  = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_LENGTH   = 400

# ── Placeholder returned when ESMFold is unavailable ─────────────────────────
def _placeholder_pdb(sequence: str, reason: str) -> str:
    return f"""HEADER    STRUCTURE PREDICTION UNAVAILABLE
REMARK    ESMFold API unavailable: {reason}
REMARK    This is a placeholder — not real coordinates.
REMARK    Rerun when ESMFold is accessible.
ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C
TER
END
"""


def _parse_plddt(pdb_text: str) -> float | None:
    """
    Parse mean pLDDT from the B-factor column of ATOM records.
    PDB B-factor column is characters 61-66 (1-indexed).
    ESMFold stores per-residue pLDDT (0-100) there.
    Returns mean pLDDT as float, or None if unparseable.
    """
    scores = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and len(line) >= 66:
            try:
                bfactor = float(line[60:66].strip())
                scores.append(bfactor)
            except ValueError:
                continue

    if not scores:
        return None

    # One B-factor per atom — take mean over CA atoms only to avoid bias
    ca_scores = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and len(line) >= 66 and " CA " in line[12:16]:
            try:
                ca_scores.append(float(line[60:66].strip()))
            except ValueError:
                continue

    values = ca_scores if ca_scores else scores
    return round(sum(values) / len(values), 2)


def _interpret_plddt(score: float | None) -> str:
    if score is None:
        return "Unavailable"
    if score >= 90:
        return "Very high confidence (pLDDT ≥ 90)"
    if score >= 70:
        return "Confident (pLDDT 70–90)"
    if score >= 50:
        return "Low confidence (pLDDT 50–70)"
    return "Very low confidence (pLDDT < 50)"


# ── Main ─────────────────────────────────────────────────────────────────────
async def run(sequence: str, **kwargs) -> dict:
    """
    Structure Prediction via ESMFold API (Meta AI).

    Returns real PDB coordinates and mean pLDDT confidence score.
    Falls back to a labelled placeholder if ESMFold is unreachable.
    """
    if not sequence:
        raise ValueError("Sequence is required for Structure Prediction.")

    # Truncate gracefully if over ESMFold limit
    if len(sequence) > MAX_LENGTH:
        print(f"  [ESMFold] Sequence length {len(sequence)} exceeds {MAX_LENGTH} AA limit — truncating.")
        sequence = sequence[:MAX_LENGTH]

    pdb_text         = None
    confidence_score = None
    source           = "ESMFold (Meta AI)"

    try:
        async with httpx.AsyncClient(timeout=120, follow_redirects=True) as client:
            resp = await client.post(
                ESMFOLD_URL,
                content=sequence.encode(),
                headers={"Content-Type": "text/plain"},
            )
            resp.raise_for_status()
            pdb_text = resp.text

        # Sanity check — real PDB has ATOM records
        if "ATOM" not in pdb_text:
            raise ValueError("ESMFold returned response with no ATOM records")

        confidence_score = _parse_plddt(pdb_text)
        print(f"  [ESMFold] Success — pLDDT={confidence_score}")

    except Exception as e:
        print(f"  [ESMFold] Failed ({e}) — returning placeholder PDB")
        pdb_text         = _placeholder_pdb(sequence, str(e))
        confidence_score = None
        source           = "Placeholder (ESMFold unavailable)"

    return {
        "predicted_pdb":           pdb_text,
        "confidence_score":        confidence_score,
        "c_score_interpretation":  _interpret_plddt(confidence_score),
        "structure_source":        source,
    }