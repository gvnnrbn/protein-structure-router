import io
import warnings
from Bio.PDB.PDBParser import PDBParser
from typing import List, Dict
import re
from unittest import result
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO

def parse_fasta_to_sequence(fasta_content: str) -> str:
    """
    Parses a FASTA string and extracts the first amino acid sequence found.
    Since the sequence will be used to search for a matching PDB structure, grabbing the first record is sufficient.
    """
    lines = fasta_content.splitlines()
    sequence_lines = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        if line.startswith(">"):
            if sequence_lines:
                break
            continue
            
        sequence_lines.append(line)
        
    return "".join(sequence_lines)

def convert_mmcif_to_pdb(mmcif_content: str) -> str:
    try:
        mmcif_file_handle = io.StringIO(mmcif_content)
        
        #  parse the mmCIF structure
        parser = MMCIFParser(QUIET=True) 
        structure = parser.get_structure("temp_structure", mmcif_file_handle)
        
        # Write PDB format to a string buffer
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        
        pdb_file_handle = io.StringIO()
        pdb_io.save(pdb_file_handle)
        
        pdb_text = pdb_file_handle.getvalue()

        pdb_id = None
        match = re.search(r'_entry\.id\s+([A-Za-z0-9_]+)', mmcif_content)
        if match:
            pdb_id = match.group(1).upper()
        elif "idcode" in structure.header and structure.header["idcode"]:
            pdb_id = structure.header["idcode"].upper()

        if pdb_id:
            header_line = f"HEADER    PROTEIN STRUCTURE                       {pdb_id}\n"
            pdb_text = header_line + pdb_text
        
        mmcif_file_handle.close()
        pdb_file_handle.close()
        print("[CONVERTER] mmCIF to PDB conversion successful:",pdb_id if pdb_id else "Unknown PDB ID")
        return {"status": "success", "data": pdb_text, "message": None, "protein_id":pdb_id if pdb_id else "Unknown", "id_type": "pdb"}

    except Exception as e:
        err_msg = f"Error converting mmCIF to PDB: {e}"
        print(f"[CONVERTER]: {err_msg}")
        return {"status": "error", "data": None, "message": err_msg,"protein_id":pdb_id if pdb_id else "Unknown", "id_type": "pdb"}

def extract_pdb_metadata_and_chains(pdb_content: str) -> dict:
    """
    Single-pass parser for PDB content.
    Extracts global metadata (protein_id, id_type) and all available chains
    with their respective sequences and lengths.
    id_type will strictly be "pdb", "uniprot", or "unknown".
    """
    protein_id = None
    id_type = "unknown"
    
    seqres_data = {}
    atom_data = {}
    last_res_num = {}
    
    three_to_one = {
        'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
        'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I',
        'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
        'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'
    }
    
    lines = pdb_content.splitlines()
    
    for line in lines:
        # 1. Extract PDB ID from HEADER (Highest Priority)
        if line.startswith("HEADER"):
            match = re.search(r'([A-Za-z0-9]{4})$', line.strip())
            if match:
                protein_id = match.group(1)
                id_type = "pdb"
                
        # 2. Extract UniProt ID from TITLE (AlphaFold specific, fallback)
        elif line.startswith("TITLE ") and id_type != "pdb":
            if "ALPHAFOLD" in line.upper():
                # Example: TITLE     ALPHAFOLD MONOMER v2.0 PREDICTION FOR O15294
                match = re.search(r'FOR\s+([A-Za-z0-9]{6,10})', line.upper())
                if match:
                    protein_id = match.group(1)
                    id_type = "uniprot"
                    
        # 3. Extract UniProt ID from DBREF (Fallback)
        elif line.startswith("DBREF") and id_type != "pdb":
            parts = line.split()
            if len(parts) > 6:
                db_name = parts[5]
                db_accession = parts[6]
                if db_name in ["UNP", "SWS"] or re.match(r'^[O,P,Q][0-9][A-Z0-9]{3}[0-9]|[A-N,R-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', db_accession):
                    protein_id = db_accession
                    id_type = "uniprot"
                    
        # Extract sequence from SEQRES (Preferred)
        elif line.startswith("SEQRES"):
            chain_id = line[11]
            if chain_id == ' ': 
                continue
                
            residues = line[19:].split()
            if chain_id not in seqres_data:
                seqres_data[chain_id] = []
                
            for res in residues:
                seqres_data[chain_id].append(three_to_one.get(res.upper(), 'X'))
                
        # Extract sequence from ATOM as a fallback
        elif line.startswith("ATOM  ") or line.startswith("HETATM"):
            # Only use CA atoms to avoid counting the same residue multiple times
            if line[12:16].strip() == "CA":
                chain_id = line[21]
                if chain_id == ' ': 
                    continue
                    
                res_name = line[17:20].strip()
                res_num = line[22:26].strip()
                
                if chain_id not in atom_data:
                    atom_data[chain_id] = []
                    last_res_num[chain_id] = None
                    
                if res_num != last_res_num[chain_id]:
                    atom_data[chain_id].append(three_to_one.get(res_name.upper(), 'X'))
                    last_res_num[chain_id] = res_num

    # Build the final chains array
    chains = []
    all_chain_ids = set(seqres_data.keys()).union(set(atom_data.keys()))
    
    for cid in sorted(all_chain_ids):
        # Prioritize SEQRES over ATOM as it contains the full intended sequence
        if cid in seqres_data and len(seqres_data[cid]) > 0:
            seq_str = "".join(seqres_data[cid])
        elif cid in atom_data and len(atom_data[cid]) > 0:
            seq_str = "".join(atom_data[cid])
        else:
            seq_str = ""
            
        if seq_str:
            chains.append({
                "chain_id": cid,
                "sequence": seq_str,
                "length": len(seq_str)
            })
            
    return {
        "protein_id": protein_id,
        "id_type": id_type,
        "chains": chains
    }