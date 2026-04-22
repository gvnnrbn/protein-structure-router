from src.validators import (
    is_valid_pdb_id,
    is_valid_alphafold_id,
    is_valid_uniprot_accession,
    is_amino_acid_sequence,
    is_fasta_format,
    is_pdb_format,
    is_mmcif_format
)
from src.fetchers import (
    search_rcsb_by_pdb_id,
    fetch_alphafold_model,
    search_rcsb_by_sequence
)
from src.converters import parse_fasta_to_sequence, convert_mmcif_to_pdb

def structure_router(query_type: str, text_query: str = None, file_content: str = None):
    if query_type == "file" and file_content:
        file_content = file_content.strip()
        
        # PDB file: pass it straight through
        if is_pdb_format(file_content):
            print("[ROUTER] PDB file detected. Passing through.")
            return file_content
            
        # mmCIF file: convert it to PDB format
        if is_mmcif_format(file_content):
            print("[ROUTER] mmCIF file detected. Converting to PDB...")
            return convert_mmcif_to_pdb(file_content)
            
        # FASTA file: extract sequence and search for pdb
        if is_fasta_format(file_content):
            print("[ROUTER] FASTA file detected. Extracting sequence...")
            raw_sequence = parse_fasta_to_sequence(file_content)
            if raw_sequence and is_amino_acid_sequence(raw_sequence):
                return search_rcsb_by_sequence(raw_sequence)
            else:
                print("[ROUTER] Error: FASTA contained invalid amino acids or was too short.")
                return None
                
        print("[ROUTER] Error: Unrecognized file format.")
        return None
    
    elif query_type in ["id", "sequence"] and text_query:
        text_query = text_query.strip()
        # ID
        if is_valid_pdb_id(text_query):
            return search_rcsb_by_pdb_id(text_query)
        
        if is_valid_alphafold_id(text_query):
            return fetch_alphafold_model(specific_af_id=text_query)
        
        if is_valid_uniprot_accession(text_query):
            return fetch_alphafold_model(uniprot_id=text_query)
        
        # SEQUENCE
        if is_amino_acid_sequence(text_query):
            return search_rcsb_by_sequence(text_query)

        print("[ROUTER] Error: Unrecognized text format.")
        return None
        
    print("[ROUTER] Error: Invalid request schema.")
    return None