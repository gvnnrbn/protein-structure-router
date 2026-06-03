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
        
        # mmCIF file: convert it to PDB format
        if is_mmcif_format(file_content):
            print("[ROUTER] mmCIF file detected. Converting to PDB...")
            result = convert_mmcif_to_pdb(file_content)
            result["format"] = "mmCIF"
            result["status"] = "success"
            return result
        
        # PDB file: pass it straight through
        if is_pdb_format(file_content):
            print("[ROUTER] PDB file detected. Passing through.")
            return {
                "status": "success", 
                "data": file_content, 
                "message": None, 
                "format": "PDB"
            }
            
        
            
        # FASTA file: extract sequence and search for pdb
        if is_fasta_format(file_content):
            print("[ROUTER] FASTA file detected. Extracting sequence...")
            raw_sequence = parse_fasta_to_sequence(file_content)
            if raw_sequence and is_amino_acid_sequence(raw_sequence):
                result = search_rcsb_by_sequence(raw_sequence)
                result["format"] = "FASTA"
                result["status"] = "success"
                return result
            else:
                err_msg = "FASTA contained invalid amino acids or was too short"
                print(f"[ROUTER] Error: {err_msg}")
                return {"status": "error", "data": None, "message": err_msg, "format": "FASTA"}
        err_msg = "Unrecognized file format."        
        print(f"[ROUTER] Error: {err_msg}")
        return {"status": "error", "data": None, "message": err_msg, "format": "Unknown"}
    
    elif query_type == "text" and text_query:
        text_query = text_query.strip()
        # ID
        if is_valid_pdb_id(text_query):
            result = search_rcsb_by_pdb_id(text_query)
            result["format"] = "PDB ID"
            result["status"] = "success"
            return result
        
        if is_valid_alphafold_id(text_query):
            result = fetch_alphafold_model(specific_af_id=text_query)
            result["format"] = "AlphaFold ID"
            result["status"] = "success"
            return result
        
        if is_valid_uniprot_accession(text_query):
            result = fetch_alphafold_model(uniprot_id=text_query)
            result["format"] = "Uniprot Accession"
            return result
        
        # SEQUENCE
        if is_amino_acid_sequence(text_query):
            result = search_rcsb_by_sequence(text_query)
            result["format"] = "Sequence"
            return result
        
        err_msg = "Unrecognized text format."
        print(f"[ROUTER] Error: {err_msg}")
        return {"status": "error", "data": None, "message": err_msg, "format": "Unknown"}

    err_msg = "Invalid request schema."    
    print(f"[ROUTER] Error: {err_msg}")
    return {"status": "error", "data": None, "message": err_msg, "format": "Unknown"}