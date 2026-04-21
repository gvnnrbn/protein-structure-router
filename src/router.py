from src.validators import (
    is_valid_pdb_id,
    is_valid_alphafold_id,
    is_valid_uniprot_accession,
    is_amino_acid_sequence,
    is_fasta_format 
)
from src.fetchers import (
    search_rcsb_by_pdb_id,
    search_alphafold_by_uniprot_id,
    search_rcsb_by_sequence
)
from src.converters import parse_fasta_to_sequence 
def structure_router(user_input: str):
    user_input = user_input.strip()
    
    if is_valid_pdb_id(user_input):
        return search_rcsb_by_pdb_id(user_input)
    
    if is_valid_alphafold_id(user_input):
        uniprot_id = user_input.split('-')[1]
        return search_alphafold_by_uniprot_id(uniprot_id)
    
    if is_valid_uniprot_accession(user_input):
        return search_alphafold_by_uniprot_id(user_input)
    
    if is_fasta_format(user_input):
        print("[ROUTER] FASTA format detected. Extracting sequence...")
        raw_sequence = parse_fasta_to_sequence(user_input)
        if raw_sequence and is_amino_acid_sequence(raw_sequence):
            return search_rcsb_by_sequence(raw_sequence)
        else:
            print("[ROUTER] Error: FASTA contained invalid amino acids or was too short.")
            return None

    # 5. Check for Raw Sequence (Fallback)
    if is_amino_acid_sequence(user_input):
        return search_rcsb_by_sequence(user_input)
        
    print("[ROUTER] Unrecognized input format.")
    return None