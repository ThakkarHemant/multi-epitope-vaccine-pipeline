from Bio.SeqUtils.ProtParam import ProteinAnalysis

async def run(sequence: str, **kwargs) -> dict:
    """
    Stage 2: Physicochemical Analysis
    Analyzes the full protein sequence for stability, MW, pI, etc.
    
    Decision Gate:
    If instability index > 40 -> mark as unstable.
    """
    if not sequence:
        raise ValueError("Sequence is required for Physicochemical Analysis.")
        
    # Clean sequence just in case (e.g. removing whitespace, unknown amino acids if needed, 
    # but BioPython usually handles standard standard AAs. We'll replace tricky ones if any,
    # but usually ProtParam expects standard 20 amino acids + some ambiguities).
    clean_seq = sequence.replace(" ", "").replace("\n", "").upper()
    
    # We strip 'X' or 'Z' or 'B' if they cause issues with ProtParam, or handle them.
    # ProtParam can raise errors on non-standard aminos.
    analysis = ProteinAnalysis(clean_seq)
    
    molecular_weight = round(analysis.molecular_weight(), 2)
    pI = round(analysis.isoelectric_point(), 2)
    instability_index = round(analysis.instability_index(), 2)
    gravy = round(analysis.gravy(), 3)
    
    # Aliphatic index requires calculation of mole fraction of Ala, Val, Ile, Leu
    # Formula: aliphatic_index = X(Ala) + a*X(Val) + b*(X(Ile) + X(Leu))
    # where a = 2.9, b = 3.9
    aa_percents = analysis.get_amino_acids_percent()
    aliphatic_index = round(
        (aa_percents.get('A', 0) * 100) + 
        2.9 * (aa_percents.get('V', 0) * 100) + 
        3.9 * ((aa_percents.get('I', 0) + aa_percents.get('L', 0)) * 100), 
        2
    )
    
    is_stable = instability_index <= 40
    
    return {
        "molecular_weight": molecular_weight,
        "theoretical_pi": pI,
        "instability_index": instability_index,
        "aliphatic_index": aliphatic_index,
        "gravy": gravy,
        "is_stable": is_stable,
        "gate_passed": is_stable
    }
