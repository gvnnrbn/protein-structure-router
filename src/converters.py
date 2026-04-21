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
        
        return pdb_text

    except Exception as e:
        print(f"[CONVERTER] Error converting mmCIF to PDB: {e}")
        return None

def clean_pdb_text(dirty_text):
    """
    Takes a dirty PDB string (with literal \\n and escape \\) 
    and converts it into a strict, valid PDB format.
    """
    # 1. Replace literal '\n' (text) with actual newlines
    if '\\n' in dirty_text:
        dirty_text = dirty_text.replace('\\n', '\n')
        
    clean_lines = []
    
    # Valid PDB keywords (Record Types) to keep: ATOM, HETATM, TER, and END.
    valid_keywords = (
        'HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 'EXPDTA', 
        'AUTHOR', 'REVDAT', 'JRNL', 'REMARK', 'DBREF', 'SEQRES', 
        'FORMUL', 'HELIX', 'SHEET', 'CRYST1', 'ORIGX', 'SCALE', 
        'ATOM', 'HETATM', 'TER', 'END', 'MASTER'
    )
    
    for line in dirty_text.split('\n'):
        # 2. Strip backslashes '\' and spaces from the START of the line.
        trimmed_line = re.sub(r'^[\s\\]+', '', line)
        
        # 3. Filter: If the line starts with a valid keyword, we keep it
        if trimmed_line.startswith(valid_keywords):
            clean_lines.append(trimmed_line)
            
    # Join everything with clean newlines
    return '\n'.join(clean_lines) + '\n'
