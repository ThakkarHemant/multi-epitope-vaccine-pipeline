async def run(
    safe_mhc1_epitopes: list,
    safe_mhc2_epitopes: list,
    bcell_epitopes: list,
    **kwargs,
) -> dict:

    RS09    = "APPHALS"
    HP91    = "KDPNAPKRPPSAFFLFCSE"
    PADRE   = "AAATLWKAAKFVA"

    EAAAK   = "EAAAK"
    GPGPG   = "GPGPG"
    AAY     = "AAY"
    KK      = "KK"

    # Select epitopes
    mhc2_peps = [e["peptide"] for e in safe_mhc2_epitopes[:2]]
    mhc1_peps = [e["peptide"] for e in safe_mhc1_epitopes[:4]]
    bcell_peps = [e["peptide"] for e in bcell_epitopes[:1]]

    parts = [RS09, EAAAK, PADRE]

    # MHC-II
    if mhc2_peps:
        parts.append(GPGPG.join(mhc2_peps))

    # MHC-I
    if mhc1_peps:
        parts.append(AAY.join(mhc1_peps))

    # B-cell
    if bcell_peps:
        parts.append(KK.join(bcell_peps))

    parts.append(EAAAK)
    parts.append(HP91)

    construct_sequence = "".join(parts)

    return {
        "adjuvant_n_terminal": RS09,
        "adjuvant_c_terminal": HP91,
        "padre_sequence": PADRE,
        "mhc2_epitopes_used": mhc2_peps,
        "mhc1_epitopes_used": mhc1_peps,
        "bcell_epitopes_used": bcell_peps,
        "final_construct_sequence": construct_sequence,
        "construct_length_aa": len(construct_sequence),
    }