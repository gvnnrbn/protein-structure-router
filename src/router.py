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
from src.converters import (
    parse_fasta_to_sequence, 
    convert_mmcif_to_pdb, 
    extract_pdb_metadata_and_chains, 
)

def structure_router(query_type: str, text_query: str = None, file_content: str = None) -> dict:
    result = None
    
    # ---------------------------------------------------------
    # 1. FETCHING AND VALIDATION 
    # ---------------------------------------------------------
    
    if query_type == "file" and file_content:
        print("[ROUTER] File upload detected. Processing...")
        file_content = file_content.strip()
        
        # mmCIF file: convert it to PDB format
        if is_mmcif_format(file_content):
            print("[ROUTER] mmCIF file detected. Converting to PDB...")
            result = convert_mmcif_to_pdb(file_content)
            result["format"] = "mmCIF"
            result["status"] = "success"
            
        # PDB file: pass it straight through
        elif is_pdb_format(file_content):
            print("[ROUTER] PDB file detected. Passing through.")
            result = {
                "status": "success", 
                "data": file_content, 
                "message": None, 
                "format": "PDB"
            }
            
        # FASTA file: extract sequence and search for pdb
        elif is_fasta_format(file_content):
            print("[ROUTER] FASTA file detected. Extracting sequence...")
            raw_sequence = parse_fasta_to_sequence(file_content) 
            if raw_sequence and is_amino_acid_sequence(raw_sequence):
                result = search_rcsb_by_sequence(raw_sequence)
                result["format"] = "FASTA"
            else:
                err_msg = "FASTA contained invalid amino acids or was too short"
                print(f"[ROUTER] Error: {err_msg}")
                return {"status": "error", "data": None, "message": err_msg, "format": "FASTA"}
                
        else:
            err_msg = "Unrecognized file format."        
            print(f"[ROUTER] Error: {err_msg}")
            return {"status": "error", "data": None, "message": err_msg, "format": "Unknown"}
    
    elif query_type == "text" and text_query:
        print("[ROUTER] Text query detected. Processing...")
        text_query = text_query.strip()
        
        # ID: PDB
        if is_valid_pdb_id(text_query):
            print("[ROUTER] PDB ID detected. Fetching structure...")
            result = search_rcsb_by_pdb_id(text_query)
            result["format"] = "PDB ID"
        
        # ID: AlphaFold Specific
        elif is_valid_alphafold_id(text_query):
            print("[ROUTER] AlphaFold ID detected. Fetching structure...")
            result = fetch_alphafold_model(specific_af_id=text_query)
            result["format"] = "AlphaFold ID"
            
        # ID: UniProt
        elif is_valid_uniprot_accession(text_query):
            print("[ROUTER] Uniprot Accession detected. Fetching AlphaFold model...")
            result = fetch_alphafold_model(uniprot_id=text_query)
            result["format"] = "Uniprot Accession"
            
        # SEQUENCE
        elif is_amino_acid_sequence(text_query):
            print("[ROUTER] Amino acid sequence detected. Searching for matching PDB structures...")
            result = search_rcsb_by_sequence(text_query)
            result["format"] = "Sequence"
            
        else:
            err_msg = "Unrecognized text format."
            print(f"[ROUTER] Error: {err_msg}")
            return {"status": "error", "data": None, "message": err_msg, "format": "Unknown"}

    else:
        err_msg = "Invalid request schema."    
        print(f"[ROUTER] Error: {err_msg}")
        return {"status": "error", "data": None, "message": err_msg, "format": "Unknown"}

    # Guard clause: if fetch failed, return the error immediately
    if result.get("status") == "error":
        return result


   # ---------------------------------------------------------
    # 2. ORQUESTRATOR
    # ---------------------------------------------------------
    
    # SCENARIO C: AlphaFold returned multiple models (multiple_choices)
    if result.get("status") == "multiple_choices":
        return {
            "status": "multiple_choices",
            "choice_type": "alphafold_models",
            "options": result.get("data", []),
            "input_format": result.get("format", "Unknown"),
            "protein_id": text_query if query_type == "text" else "unknown",
            "id_type": "uniprot",
            "pdb_found": None, 
        }
        
    # SCENARIOS A and B: We have the raw PDB string in our hands (success)
    elif result.get("status") == "success":
        raw_pdb = result.get("data")
        
        # ONE SINGLE PASS: Extract IDs, chains, sequences, and lengths
        metadata = extract_pdb_metadata_and_chains(raw_pdb)
        detected_chains = metadata.get("chains", [])
        
        # Prioritize extracted IDs. Fallback to text_query only if extraction fails.
        extracted_protein_id = metadata.get("protein_id")
        if not extracted_protein_id:
            extracted_protein_id = text_query if query_type == "text" else "unknown"
            
        extracted_id_type = metadata.get("id_type", "unknown")
        
        # SCENARIO B: Multiple chains in a single file/download
        if len(detected_chains) > 1:
            print(f"[ROUTER] Multiple chains detected ({len(detected_chains)}). Prompting user.")
            
            return {
                "status": "multiple_choices",
                "choice_type": "chains",        
                "options": detected_chains,      
                "input_format": result.get("format", "Unknown"),
                "protein_id": extracted_protein_id,
                "id_type": extracted_id_type,
                "pdb_found": raw_pdb,
            }
            
        # SCENARIO A: Single chain detected
        else:
            print("[ROUTER] Single chain detected. Ready for immediate TAPO execution.")
            
            single_chain_info = detected_chains[0] if detected_chains else {
                "chain_id": "A", 
                "sequence": "", 
                "length": 0
            }
            
            return {
                "status": "success",
                "input_format": result.get("format", "Unknown"),
                "protein_id": extracted_protein_id,
                "id_type": extracted_id_type,
                "chain_id": single_chain_info["chain_id"],
                "length": single_chain_info["length"],
                "sequence": single_chain_info["sequence"],
                "pdb_found": raw_pdb,
            }

    return result