import httpx
import re


async def run(sequence: str, uniprot_id: str, **kwargs) -> dict:
    """
    Step 1: Sequence Retrieval
    Fetches the protein sequence and metadata from UniProt.
    Supports canonical IDs (e.g. Q5JQC9) and isoform IDs (e.g. Q5JQC9-1, Q5JQC9-2).
    """
    if not uniprot_id:
        raise ValueError("UniProt ID is required for Sequence Retrieval.")

    uid = uniprot_id.strip().upper()
    is_isoform = bool(re.match(r'^[A-Z0-9]+-\d+$', uid))

    async with httpx.AsyncClient(timeout=30.0) as client:
        if is_isoform:
            # Isoform IDs: use FASTA endpoint to get variant sequence
            base_id = uid.split("-")[0]
            fasta_url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
            fasta_resp = await client.get(fasta_url)
            if fasta_resp.status_code != 200:
                raise Exception(f"Failed to fetch isoform FASTA from UniProt for {uid}. "
                                f"Status: {fasta_resp.status_code}")
            fasta_text = fasta_resp.text
            lines = fasta_text.strip().splitlines()
            header = lines[0] if lines else ""
            protein_seq = "".join(lines[1:])

            # Parse protein name from FASTA header e.g. >sp|Q5JQC9-1|... Protein name OS=...
            protein_name = "Unknown Protein"
            m = re.search(r'\|[^|]+\|\S+\s+(.+?)\s+OS=', header)
            if m:
                protein_name = m.group(1)

            # Get organism from base entry JSON
            json_url = f"https://rest.uniprot.org/uniprotkb/{base_id}.json"
            json_resp = await client.get(json_url)
            organism = "Unknown Organism"
            if json_resp.status_code == 200:
                organism = json_resp.json().get("organism", {}).get("scientificName", "Unknown Organism")

        else:
            # Canonical ID: use JSON API
            url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
            response = await client.get(url)
            if response.status_code != 200:
                raise Exception(f"Failed to fetch from UniProt for {uid}. "
                                f"Status: {response.status_code}")
            data = response.json()
            protein_seq = data.get("sequence", {}).get("value", "")
            if not protein_seq:
                raise Exception("Sequence not found in UniProt response.")
            protein_names = data.get("proteinDescription", {}).get("recommendedName", {})
            protein_name = protein_names.get("fullName", {}).get("value", "Unknown Protein")
            organism = data.get("organism", {}).get("scientificName", "Unknown Organism")

    return {
        "uniprot_id": uid,
        "isoform": is_isoform,
        "sequence": protein_seq,
        "length": len(protein_seq),
        "protein_name": protein_name,
        "organism": organism,
    }
