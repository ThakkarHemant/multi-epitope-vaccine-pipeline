import hashlib

# HLA coverage bases per region
COVERAGE_REGIONS = {
    "North America": 0.89,
    "Europe": 0.87,
    "East Asia": 0.82,
    "South Asia": 0.85,
    "Africa": 0.78,
    "World": 0.84,
}


async def run(bcell_epitopes: list, mhc1_epitopes: list, mhc2_epitopes: list, **kwargs) -> dict:
    """
    Step 6: Population Coverage Analysis (IEDB stub)
    Input: Safe epitopes from Step 5 (B-cell, MHC-I, MHC-II).
    Output: Coverage percentages per region and shortlisted alleles.
    """
    bcells = [e["peptide"] for e in bcell_epitopes]
    mhc1_alleles = list({e["hla_allele"] for e in mhc1_epitopes})
    mhc2_alleles = list({e["hla_allele"] for e in mhc2_epitopes})

    # Deterministic regional coverage variation based on inputs
    seed_str = "".join(bcells + mhc1_alleles + mhc2_alleles)
    seed = int(hashlib.md5(seed_str.encode()).hexdigest(), 16)
    variation = (seed % 10) / 100.0  # ±0–0.09
    
    # B-cell coverage is usually broader globally, simulated here:
    bcell_variation = ((seed // 10) % 5) / 100.0
    
    bcell_coverage = {r: round(min(0.99, base + 0.05 + bcell_variation), 4) for r, base in COVERAGE_REGIONS.items()}
    class1_coverage = {r: round(min(0.99, base + variation), 4) for r, base in COVERAGE_REGIONS.items()}
    class2_coverage = {r: round(min(0.99, base - 0.04 + variation), 4) for r, base in COVERAGE_REGIONS.items()}

    return {
        "bcell_peptides_used": bcells,
        "mhc1_alleles_used": mhc1_alleles,
        "mhc2_alleles_used": mhc2_alleles,
        "bcell_population_coverage": bcell_coverage,
        "class1_population_coverage": class1_coverage,
        "class2_population_coverage": class2_coverage,
        "world_bcell_coverage_pct": round(bcell_coverage["World"] * 100, 2),
        "world_class1_coverage_pct": round(class1_coverage["World"] * 100, 2),
        "world_class2_coverage_pct": round(class2_coverage["World"] * 100, 2),
    }
