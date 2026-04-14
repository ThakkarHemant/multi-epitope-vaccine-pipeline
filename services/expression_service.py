import asyncio
import hashlib

async def run(sequence: str, host: str = "E. coli", **kwargs) -> dict:
    """
    Stage 10: Codon Optimization
    Simulates JCat codon optimization for a specific expression host.
    
    Returns optimized DNA sequence, CAI score, and GC content.
    """
    if not sequence:
        raise ValueError("Protein sequence is required for Codon Optimization.")
        
    await asyncio.sleep(0.5)
    
    # Mapping amino acids to highly optimized codons for E. coli (simplistic stub)
    codon_map_ecoli = {
        'A': 'GCG', 'R': 'CGC', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
        'Q': 'CAG', 'E': 'GAA', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT',
        'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
        'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
        '*': 'TAA'
    }
    
    clean_seq = sequence.replace(" ", "").replace("\n", "").upper()
    
    optimized_dna = []
    for aa in clean_seq:
        # Default to a generic codon if unusual AA
        codon = codon_map_ecoli.get(aa, 'NNN')
        optimized_dna.append(codon)
        
    dna_string = "".join(optimized_dna)
    
    # Calculate actual GC content of our stub sequence
    gc_count = dna_string.count('G') + dna_string.count('C')
    gc_content = round((gc_count / len(dna_string)) * 100, 2) if dna_string else 0.0
    
    # Codon Adaptation Index (CAI) ideally > 0.8
    # Based on our static map to optimal codons, we simulate a CAI around 0.95
    hash_val = int(hashlib.md5(sequence.encode("utf-8")).hexdigest(), 16)
    cai_score = round(0.90 + ((hash_val % 10) / 100.0), 2)
    
    return {
        "host_organism": host,
        "optimized_dna": dna_string,
        "cai_score": cai_score,
        "gc_content": gc_content
    }
