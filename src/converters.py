import io
import re
from unittest import result
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBParser
from src.validators import (
    is_valid_pdb_id,
    is_valid_uniprot_accession,
    is_pdb_format,
)


def parse_fasta_to_sequence(fasta_content: str, chain_id: str) -> str:
    lines = fasta_content.strip().split('\n')
    
    sequences = {}
    current_header = ""
    
    # Group sequence blocks by their FASTA headers
    for line in lines:
        clean_line = line.strip()
        if not clean_line:
            continue
            
        if clean_line.startswith('>'):
            current_header = clean_line
            sequences[current_header] = []
        elif current_header:
            sequences[current_header].append(clean_line)
        else:
            # Handles raw sequence files without any '>' headers
            current_header = "raw_sequence"
            sequences[current_header] = [clean_line]
            
    target_sequence_blocks = []
    
    # Find the sequence that matches the requested chain_id
    for header, blocks in sequences.items():
        if header == "raw_sequence":
            target_sequence_blocks = blocks
            break
            
        # Matches standard RCSB PDB format: "|Chains A, C[auth D]|"
        chain_match = re.search(r'Chains?\s+([^|]+)', header, re.IGNORECASE)
        
        if chain_match:
            valid_chains = re.findall(r'\b[A-Za-z0-9]+\b', chain_match.group(1))
            if chain_id in valid_chains:
                target_sequence_blocks = blocks
                break
        else:
            # Fallback for generic FASTA files where the ID is directly in the header
            if f">{chain_id}" in header or f"|{chain_id}|" in header:
                target_sequence_blocks = blocks
                break
                
    # Fallback to the first available sequence if no chain strictly matches
    if not target_sequence_blocks and sequences:
        target_sequence_blocks = list(sequences.values())[0]

    raw_sequence = "".join(target_sequence_blocks).upper()
    
    if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', raw_sequence):
        print(f"[CONVERTER] Warning: FASTA sequence for chain {chain_id} contains non-standard amino acids or gaps.")
        
    return raw_sequence

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

def clean_pdb_text(dirty_text: str) -> str:
    if '\\n' in dirty_text:
        dirty_text = dirty_text.replace('\\n', '\n')
        
    clean_lines = []
    is_capturing = False
    
    for line in dirty_text.split('\n'):
        trimmed_line = re.sub(r'^[\s\\]+', '', line)
        clean_lines.append(trimmed_line)
                
    return '\n'.join(clean_lines) + '\n' if clean_lines else ''

def get_chains_from_pdb(pdb_string: str) -> list:
    """
    Extracts a list of available chains (e.g., ['A', 'B']) from a PDB string using BioPython.
    """
    # QUIET=True prevents warnings for imperfect or non-standard PDB formats
    parser = PDBParser(QUIET=True) 
    structure = parser.get_structure("temp_structure", io.StringIO(pdb_string))
    
    chains = []
    for model in structure:
        for chain in model:
            if chain.id not in chains:
                chains.append(chain.id)
        break 
    return chains



def extract_pdb_metadata(result: dict, chain_id: str) -> dict:
    """
    Extracts protein ID, ID type (PDB or UniProt), sequence, and length from the PDB content.
    """
    if not result or "data" not in result or not result["data"]:
        return result
        
    pdb_content = result["data"]
    
    # Validate that it is actually a PDB format before parsing
    if not is_pdb_format(pdb_content):
        return result
        
    # 1. Extract ID and its type (PDB or UniProt)
    protein_id = None
    id_type = "unknown"
    
    lines = pdb_content.splitlines()
    
    # Iterate over the first 100 lines where metadata is usually located
    for line in lines[:100]: 
        if line.startswith("HEADER"):
            parts = line.split()
            if parts:
                candidate = parts[-1]
                if is_valid_pdb_id(candidate):
                    protein_id = candidate
                    id_type = "pdb"
                    break
        elif line.startswith("DBREF"):
            # DBREF usually has the PDB ID at the start and the UniProt reference in following columns
            parts = line.split()
            for part in parts:
                if is_valid_uniprot_accession(part):
                    protein_id = part
                    id_type = "uniprot"
                    break
        elif line.startswith("TITLE") or line.startswith("COMPND"):
            # Fallback for AlphaFold files where the UniProt ID might be in Title/Compnd
            for word in line.split():
                # Remove attached punctuation marks just in case
                clean_word = ''.join(c for c in word if c.isalnum() or c == '-')
                if is_valid_uniprot_accession(clean_word):
                    protein_id = clean_word
                    id_type = "uniprot"
                    break
                elif is_valid_pdb_id(clean_word):
                    protein_id = clean_word
                    id_type = "pdb"
                    break
                    
        if protein_id:
            break

    # 2. Extract sequence based on the selected chain (chain_id)
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    sequence = []
    
    # Try with SEQRES first (the official and recommended method)
    for line in lines:
        if line.startswith("SEQRES"):
            parts = line.split()
            # Column 3 in SEQRES contains the chain identifier
            if len(parts) > 2 and parts[2] == chain_id:
                for res in parts[4:]:
                    sequence.append(three_to_one.get(res, 'X'))
                    
    # If SEQRES is missing, iterate through ATOMs as a fallback
    if not sequence:
        last_res_seq = None
        for line in lines:
            if line.startswith("ATOM"):
                # In PDB format, the chain identifier is located at column 22 (index 21)
                chain = line[21]
                if chain == chain_id:
                    res_name = line[17:20].strip()
                    res_seq = line[22:26].strip()
                    
                    # Avoid adding the same amino acid multiple times for each of its atoms
                    if res_seq != last_res_seq:
                        sequence.append(three_to_one.get(res_name, 'X'))
                        last_res_seq = res_seq
                        
    seq_str = "".join(sequence)
    
    # 3. Assign the new values to the result
    result["protein_id"] = protein_id
    result["id_type"] = id_type
    result["sequence"] = seq_str
    result["length"] = len(seq_str)
    
    return result