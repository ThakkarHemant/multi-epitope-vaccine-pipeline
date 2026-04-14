from Bio.SeqUtils.ProtParam import ProteinAnalysis


async def run(sequence: str, **kwargs) -> dict:
    """
    Step 3: Physicochemical Analysis (ProtParam)
    Returns MW, pI, instability index, aliphatic index, GRAVY.
    Informational only — no gate.
    """
    if not sequence:
        raise ValueError("Sequence is required for Physicochemical Analysis.")

    clean_seq = sequence.replace(" ", "").replace("\n", "").upper()
    analysis = ProteinAnalysis(clean_seq)

    molecular_weight = round(analysis.molecular_weight(), 2)
    pI = round(analysis.isoelectric_point(), 2)
    instability_index = round(analysis.instability_index(), 2)
    gravy = round(analysis.gravy(), 3)

    aa_percents = analysis.get_amino_acids_percent()
    aliphatic_index = round(
        (aa_percents.get("A", 0) * 100)
        + 2.9 * (aa_percents.get("V", 0) * 100)
        + 3.9 * ((aa_percents.get("I", 0) + aa_percents.get("L", 0)) * 100),
        2,
    )

    return {
        "molecular_weight": molecular_weight,
        "theoretical_pi": pI,
        "instability_index": instability_index,
        "aliphatic_index": aliphatic_index,
        "gravy": gravy,
    }
