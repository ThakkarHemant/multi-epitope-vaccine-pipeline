"""
pepfold_service.py — Step 8: Peptide 3D Structure Prediction

Primary:  ESMFold API (Meta AI) — replaces PEP-FOLD3 REST API which no longer
          exists (moved to Mobyle2 portal, no programmatic access).
          https://api.esmatlas.com/foldSequence/v1/pdb/
          - Free, no API key, accepts 5–400 AA sequences
          - Returns real PDB with per-residue pLDDT in B-factor column
          - Validated on short peptides (Lin et al. 2023, Science 379:1123)

Fallback: If ESMFold is unavailable, return a clearly-labelled placeholder PDB
          so the pipeline never crashes. confidence_score = None on fallback.

Note:  PEP-FOLD3/4 REST API was deprecated — both now require browser-based
       Mobyle2 portal submission with no documented programmatic interface.
       ESMFold is the best available free replacement for short peptide folding.

References:
  Lin Z et al. (2023) Science 379(6637):1123-30.  [ESMFold]
  Lamiable A et al. (2016) Nucleic Acids Res 44(W1):W449-54.  [PEP-FOLD3]
"""

import asyncio
import logging
import re
import httpx

logger = logging.getLogger(__name__)

ESMFOLD_URL   = "https://api.esmatlas.com/foldSequence/v1/pdb/"
VALID_AA      = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LEN, MAX_LEN = 5, 400


def _is_valid_sequence(seq: str) -> bool:
    return MIN_LEN <= len(seq) <= MAX_LEN and all(c in VALID_AA for c in seq.upper())


def _parse_plddt(pdb_text: str) -> float | None:
    """Parse mean pLDDT from CA atom B-factor column."""
    scores = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and len(line) >= 66 and line[12:16].strip() == "CA":
            try:
                scores.append(float(line[60:66].strip()))
            except ValueError:
                continue
    if not scores:
        return None
    return round(sum(scores) / len(scores), 2)


def _interpret_plddt(score: float | None) -> str:
    if score is None:       return "Unavailable"
    if score >= 90:         return "Very high confidence (pLDDT ≥ 90)"
    if score >= 70:         return "Confident (pLDDT 70–90)"
    if score >= 50:         return "Low confidence (pLDDT 50–70)"
    return "Very low confidence (pLDDT < 50)"


def _placeholder_pdb(peptide: str, reason: str) -> str:
    return (
        f"HEADER    STRUCTURE UNAVAILABLE\n"
        f"REMARK    ESMFold failed: {reason}\n"
        f"REMARK    Peptide: {peptide}\n"
        f"ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        f"TER\nEND\n"
    )


async def _fold_peptide(
    client: httpx.AsyncClient, peptide: str
) -> tuple[str, float | None, str]:
    """
    Fold a single peptide via ESMFold.
    Returns (pdb_text, plddt_score, source_label).
    """
    try:
        resp = await client.post(
            ESMFOLD_URL,
            content=peptide.encode(),
            headers={"Content-Type": "text/plain"},
            timeout=120,
        )
        resp.raise_for_status()
        pdb = resp.text

        if "ATOM" not in pdb:
            raise ValueError("No ATOM records in ESMFold response")

        plddt = _parse_plddt(pdb)
        logger.info("ESMFold: %s — pLDDT=%.1f", peptide, plddt or 0)
        return pdb, plddt, "ESMFold (Meta AI)"

    except Exception as e:
        logger.warning("ESMFold failed for %s — %s", peptide, e)
        return _placeholder_pdb(peptide, str(e)), None, f"Placeholder (ESMFold failed: {e})"


async def run(safe_mhc1_epitopes: list, **kwargs) -> dict:
    """
    Step 8: Peptide 3D Structure Prediction via ESMFold.
    Folds each safe MHC-I epitope individually.
    """
    models  = []
    errors  = []

    eligible = [
        ep for ep in safe_mhc1_epitopes
        if _is_valid_sequence(ep.get("peptide", ""))
    ]
    skipped = [
        ep["peptide"] for ep in safe_mhc1_epitopes
        if not _is_valid_sequence(ep.get("peptide", ""))
    ]

    if not eligible:
        return {
            "peptides_submitted":    0,
            "successful_models":     0,
            "pepfold_models":        [],
            "skipped_invalid_length": skipped,
            "errors":                ["No eligible peptides"],
        }

    async with httpx.AsyncClient(follow_redirects=True) as client:
        for ep in eligible:
            peptide = ep["peptide"].upper()
            hla     = ep.get("hla_allele", "unknown")

            try:
                pdb, plddt, source = await _fold_peptide(client, peptide)
                models.append({
                    "peptide":        peptide,
                    "hla_allele":     hla,
                    "pdb":            pdb,
                    "sopep_energy":   None,          # not available from ESMFold
                    "plddt_score":    plddt,
                    "confidence":     _interpret_plddt(plddt),
                    "source":         source,
                })
                if plddt is not None:
                    logger.info("ESMFold: %s modeled (pLDDT=%.1f)", peptide, plddt)
                else:
                    errors.append({"peptide": peptide, "error": "ESMFold unavailable — placeholder used"})

            except Exception as exc:
                logger.warning("Folding failed for %s — %s", peptide, exc)
                errors.append({"peptide": peptide, "error": str(exc)})

            await asyncio.sleep(1.0)   # polite delay between requests

    successful = sum(1 for m in models if m["plddt_score"] is not None)

    return {
        "peptides_submitted":    len(eligible),
        "successful_models":     successful,
        "pepfold_models":        models,
        "skipped_invalid_length": skipped,
        "errors":                errors,
    }