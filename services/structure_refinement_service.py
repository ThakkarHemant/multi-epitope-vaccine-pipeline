"""
structure_refinement_service.py — Structure Refinement

Primary:  3Drefine web server (University of Missouri)
          https://3drefine.mu.hekademeia.org/
          Bhattacharya D et al. (2016) Nucleic Acids Res 44(W1):W406-9.

          Protocol:
            1. POST multipart form to /queue/  → receive job ID
            2. Poll /queue/<job_id>/status/    → wait for FINISHED
            3. GET /queue/<job_id>/results/    → parse refined PDB + energy score

          3Drefine is async (job queue). Typical runtime < 5 min for ~200 AA.
          Returns 5 refined models; we take model 1 (best by default sort).
          Energy score: potential energy after minimization (lower = better).

Fallback: If 3Drefine is unavailable or times out, return the input PDB
          unchanged with energy_score = None and refinement_status = "Skipped".
          Pipeline never crashes.

References:
  Bhattacharya D et al. (2016) Nucleic Acids Res 44(W1):W406-9.  [3Drefine]
  Bhattacharya D & Cheng J (2013) PLoS ONE 8(7):e69648.          [i3Drefine]
"""

import asyncio
import re
import uuid
import httpx

REFINE_BASE    = "https://3drefine.mu.hekademeia.org"
SUBMIT_URL     = f"{REFINE_BASE}/queue/"
POLL_INTERVAL  = 15    # seconds between status checks
MAX_WAIT       = 600   # 10 minutes max before giving up


# ── Helpers ───────────────────────────────────────────────────────────────────

def _parse_energy(text: str) -> float | None:
    """
    Parse potential energy from 3Drefine results page or PDB REMARK lines.
    3Drefine reports energy in the results table and in REMARK records.
    Typical values: large negative floats (e.g. -15000 to -50000 kcal/mol).
    """
    # Try REMARK in PDB first
    for line in text.splitlines():
        if "REMARK" in line and "ENERGY" in line.upper():
            match = re.search(r"(-?\d+\.?\d*)", line)
            if match:
                return float(match.group(1))

    # Try HTML results table
    match = re.search(r"potential\s+energy[^-\d]*(-?\d+\.?\d*)", text, re.IGNORECASE)
    if match:
        return float(match.group(1))

    match = re.search(r"energy[^-\d]*(-\d+\.?\d*)", text, re.IGNORECASE)
    if match:
        return float(match.group(1))

    return None


def _extract_pdb(text: str) -> str | None:
    """Extract PDB content from an HTML page — look for ATOM record blocks."""
    # 3Drefine may embed PDB in <pre> tags or as plain text in a download link
    pre_match = re.search(r"<pre[^>]*>(.*?)</pre>", text, re.DOTALL | re.IGNORECASE)
    if pre_match:
        candidate = pre_match.group(1).strip()
        if "ATOM" in candidate:
            return candidate

    # If the response itself is a PDB file (direct download)
    if text.strip().startswith(("ATOM", "HEADER", "REMARK")):
        return text

    return None


# ── Main ─────────────────────────────────────────────────────────────────────

async def run(predicted_pdb: str, **kwargs) -> dict:
    """
    Structure Refinement via 3Drefine web server.

    Submits the predicted PDB, polls until job is FINISHED, then
    downloads and returns the best refined model + energy score.
    Falls back to returning the input PDB unchanged if unavailable.
    """
    if not predicted_pdb:
        raise ValueError("Predicted PDB structure is required for Refinement.")

    # Skip refinement if input is a placeholder (no real ATOM records)
    atom_lines = [l for l in predicted_pdb.splitlines() if l.startswith("ATOM")]
    if len(atom_lines) < 4:
        print("  [3Drefine] Input PDB is a placeholder — skipping refinement.")
        return {
            "refined_pdb":        predicted_pdb,
            "energy_score":       None,
            "refinement_status":  "Skipped (placeholder input)",
        }

    job_name = f"mev_{uuid.uuid4().hex[:8]}"
    refined_pdb    = None
    energy_score   = None
    status_label   = "Failed"

    try:
        async with httpx.AsyncClient(timeout=60, follow_redirects=True) as client:

            # ── Step 1: Submit job ───────────────────────────────────────────
            files = {
                "job_name":   (None, job_name),
                "structure":  ("model.pdb", predicted_pdb, "chemical/x-pdb"),
                "molprobity": (None, "on"),
                "rwplus":     (None, "on"),
            }
            resp = await client.post(SUBMIT_URL, files=files)
            resp.raise_for_status()

            # Extract job ID from redirect URL or response HTML
            job_id = None
            job_url_match = re.search(
                r"/queue/([a-zA-Z0-9_\-]+)/", resp.url.path + resp.text
            )
            if job_url_match:
                job_id = job_url_match.group(1)

            if not job_id:
                raise ValueError(f"Could not extract job ID from 3Drefine response")

            print(f"  [3Drefine] Job submitted: {job_id}")

            # ── Step 2: Poll for completion ──────────────────────────────────
            status_url  = f"{REFINE_BASE}/queue/{job_id}/status/"
            results_url = f"{REFINE_BASE}/queue/{job_id}/"
            elapsed     = 0

            while elapsed < MAX_WAIT:
                await asyncio.sleep(POLL_INTERVAL)
                elapsed += POLL_INTERVAL

                status_resp = await client.get(status_url)
                status_text = status_resp.text.strip().upper()
                print(f"  [3Drefine] Job {job_id} status: {status_text} ({elapsed}s)")

                if "FINISHED" in status_text or "COMPLETE" in status_text:
                    break
                if "ERROR" in status_text or "FAILED" in status_text:
                    raise RuntimeError(f"3Drefine job {job_id} failed: {status_text}")
            else:
                raise TimeoutError(f"3Drefine job {job_id} timed out after {MAX_WAIT}s")

            # ── Step 3: Download best refined model (model 1) ────────────────
            # Try direct PDB download first
            pdb_url  = f"{REFINE_BASE}/queue/{job_id}/model1.pdb"
            pdb_resp = await client.get(pdb_url)

            if pdb_resp.status_code == 200 and "ATOM" in pdb_resp.text:
                refined_pdb  = pdb_resp.text
                energy_score = _parse_energy(pdb_resp.text)
            else:
                # Fall back to scraping results page
                results_resp = await client.get(results_url)
                refined_pdb  = _extract_pdb(results_resp.text)
                energy_score = _parse_energy(results_resp.text)

            if not refined_pdb:
                raise ValueError("Could not extract refined PDB from 3Drefine results")

            status_label = "Successful"
            print(f"  [3Drefine] Refinement complete. Energy={energy_score}")

    except Exception as e:
        print(f"  [3Drefine] Failed ({e}) — returning input PDB unchanged")
        refined_pdb  = predicted_pdb
        energy_score = None
        status_label = f"Skipped (3Drefine unavailable: {e})"

    return {
        "refined_pdb":        refined_pdb,
        "energy_score":       energy_score,
        "refinement_status":  status_label,
    }