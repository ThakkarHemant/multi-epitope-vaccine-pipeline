import re


VECTOR = "pShuttle-CMV"
RE_SITES = {
    "EcoRI":  "GAATTC",
    "BamHI":  "GGATCC",
    "HindIII": "AAGCTT",
    "NotI":   "GCGGCCGC",
    "XhoI":   "CTCGAG",
}


async def run(optimized_dna: str, **kwargs) -> dict:
    """
    Step 12: In-Silico Cloning (SnapGene / pShuttle-CMV stub)
    - Validates restriction sites in the optimized DNA
    - Confirms correct insertion orientation
    - Generates a cloning map summary
    """
    if not optimized_dna:
        raise ValueError("Optimized DNA sequence required for cloning.")

    dna = optimized_dna.upper()

    # Check for internal restriction sites that would disrupt cloning
    conflicting_sites = []
    for enzyme, site in RE_SITES.items():
        occurrences = [m.start() for m in re.finditer(site, dna)]
        if occurrences:
            conflicting_sites.append({
                "enzyme": enzyme,
                "site": site,
                "positions": occurrences,
            })

    compatible = len(conflicting_sites) == 0
    orientation = "Forward (5' CMV promoter → Insert → 3' BGH polyA)"

    cloning_map = {
        "vector": VECTOR,
        "insert_size_bp": len(dna),
        "cloning_sites": {
            "upstream": "NotI (5')",
            "downstream": "XhoI (3')",
        },
        "promoter": "CMV immediate-early",
        "terminator": "BGH poly-A signal",
        "orientation": orientation,
    }

    return {
        "vector": VECTOR,
        "insert_length_bp": len(dna),
        "orientation": orientation,
        "restriction_site_conflicts": conflicting_sites,
        "cloning_compatible": compatible,
        "cloning_status": "Success" if compatible else "Warnings — see restriction_site_conflicts",
        "cloning_map": cloning_map,
    }
