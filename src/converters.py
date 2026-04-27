# Convert structure files into acceptable detector input
import io
import re
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO

def parse_fasta_to_sequence(fasta_content: str) -> str:
    lines = fasta_content.strip().split('\n')
    sequence_blocks = []
    
    for line in lines:
        clean_line = line.strip()
        # Ignore empty lines and header lines
        if not clean_line or clean_line.startswith('>'):
            continue
        sequence_blocks.append(clean_line)
        
    raw_sequence = "".join(sequence_blocks).upper()
    #TODO: Warning for incomplete sequencce: "X"
    if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', raw_sequence):
        print("[CONVERTER] Warning: FASTA sequence contains non-standard amino acids.")
        
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
        
        mmcif_file_handle.close()
        pdb_file_handle.close()
        
        return {"status": "success", "data": pdb_text, "message": None}

    except Exception as e:
        err_msg = f"Error converting mmCIF to PDB: {e}"
        print(f"[CONVERTER]: {err_msg}")
        return {"status": "error", "data": None, "message": err_msg}

def clean_pdb_text(dirty_text: str) -> str:
    if '\\n' in dirty_text:
        dirty_text = dirty_text.replace('\\n', '\n')
        
    clean_lines = []
    is_capturing = False
    
    for line in dirty_text.split('\n'):
        trimmed_line = re.sub(r'^[\s\\]+', '', line)
        
        if not is_capturing and trimmed_line.startswith('ATOM'):
            is_capturing = True
            
        if is_capturing:
            clean_lines.append(trimmed_line)
            if trimmed_line.startswith('TER'):
                break
                
    return '\n'.join(clean_lines) + '\n' if clean_lines else ''
