import asyncio
import hashlib


HLA_STRUCTURES = {
    "HLA-A*02:01": "PDB:1HHH",
    "HLA-B*07:02": "PDB:1BBB",
    "HLA-DRB1*03:01": "PDB:2SEB",
}


async def run(safe_mhc1_epitopes: list, safe_mhc2_epitopes: list, **kwargs) -> dict:
    """
    Step 8: Molecular Docking (HDOCK stub)
    Docks each selected epitope against its corresponding HLA structure.
    Returns docking score, binding energy, and complex identifier.
    """
    await asyncio.sleep(0.5)

    docked = []
    for ep in safe_mhc1_epitopes + safe_mhc2_epitopes:
        pep = ep["peptide"]
        allele = ep.get("hla_allele", "HLA-A*02:01")
        val = int(hashlib.md5((pep + allele).encode()).hexdigest(), 16)

        # HDOCK: docking score (higher is better, range 0–1000)
        docking_score = 300 + (val % 500)
        # Binding energy (kcal/mol — more negative = stronger binding)
        binding_energy = round(-6.0 - (val % 100) / 20.0, 2)

        docked.append({
            "peptide": pep,
            "hla_allele": allele,
            "hla_structure": HLA_STRUCTURES.get(allele, "PDB:MODELLED"),
            "docking_score": docking_score,
            "binding_energy_kcal_mol": binding_energy,
            "complex_id": f"DOCK_{pep[:4]}_{allele.replace('*','').replace(':','')}",
        })

    # Sort best to worst (highest docking score)
    docked = sorted(docked, key=lambda x: -x["docking_score"])

    return {
        "docked_complexes": docked,
        "total_docked": len(docked),
    }
