from Bio.SeqUtils.ProtParam import ProteinAnalysis

async def run(sequence: str, **kwargs) -> dict:
    """
    Stage 3: Secondary Structure Prediction
    Provides Helix %, Sheet %, Coil %.
    
    This uses BioPython's ProteinAnalysis local approximation 
    as a fast alternative to PSIPRED.
    """
    if not sequence:
        raise ValueError("Sequence is required for Secondary Structure Prediction.")
        
    clean_seq = sequence.replace(" ", "").replace("\n", "").upper()
    analysis = ProteinAnalysis(clean_seq)
    
    # Returns a tuple of (helix, turn, sheet) fractions
    ss_fractions = analysis.secondary_structure_fraction()
    
    helix_pct = round(ss_fractions[0] * 100, 2)
    sheet_pct = round(ss_fractions[2] * 100, 2)
    
    # We group Turn with Coil for simplicity or just estimate coil 
    # as 100 - (helix + sheet), since the fractions don't inherently cover everything (only 3 specific tendencies).
    # To be mathematically consistent to 100%:
    coil_pct = round(100.0 - (helix_pct + sheet_pct), 2)
    # Ensure no negative values due to floating math
    coil_pct = max(0.0, coil_pct)
    
    return {
        "helix_percent": helix_pct,
        "sheet_percent": sheet_pct,
        "coil_percent": coil_pct,
        "passed": True  # informational only
    }
