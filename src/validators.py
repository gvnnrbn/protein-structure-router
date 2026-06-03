# Identifies type of string input

import re

def is_valid_pdb_id(pdb_id: str) -> bool:
    """
    Validates a PDB ID supporting BOTH legacy and future extended formats.
    - Legacy Format: 4 alphanumeric characters, starting with 1-9 (e.g., 1LM5).
    - Extended Format (2028+): 'pdb_' prefix followed by 8 alphanumeric chars (e.g., pdb_1000axyz).
    """
    pattern = r'^([1-9][a-zA-Z0-9]{3}|pdb_[a-zA-Z0-9]{8})$'
    return bool(re.match(pattern, pdb_id, re.IGNORECASE))

def is_valid_uniprot_accession(accession: str) -> bool:
    """
    Format: 6 or 10 alphanumeric characters. Supports isoforms using -[0-9] suffix.
    """
    pattern = r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$'
    
    return bool(re.match(pattern, accession.upper()))

def is_valid_alphafold_id(af_id: str) -> bool:
    """
    Format: AF-<UniProtAC>-F<ModelVersion>
    """
    pattern = r'^AF-[A-Z0-9]{6,10}(?:-[0-9]+)?-F[0-9]+$'
    return bool(re.match(pattern, af_id))

def is_amino_acid_sequence(text: str) -> bool:
    """
    Checks if a string consists only of valid standard amino acid codes.
    """
    text = text.upper().strip()
    if len(text) < 15: 
        print("[VALIDATOR] Warning: Sequence is very short, may not be a valid protein sequence.")
        return False
    
    # Standard 20 amino acids
    aa_alphabet = set("ACDEFGHIKLMNPQRSTVWY")
    return all(char in aa_alphabet for char in text)

def is_fasta_format(text: str) -> bool:
    """
    Checks if the input string looks like a FASTA format, starting with '>' .
    """
    clean_text = text.strip()
    return clean_text.startswith(">") and "\n" in clean_text

def is_pdb_format(text: str) -> bool:
    """
    Checks if the text is a PDB file.
    """
    headers = ('HEADER', 'TITLE', 'REMARK', 'CRYST1')
    return any(text.startswith(h) for h in headers) or 'ATOM  ' in text

def is_mmcif_format(text: str) -> bool:
    """
    Checks if the text is an mmCIF file. mmCIF files must start with 'data_' and contain 'loop_'
    """
    clean_text = text.strip()
    return clean_text.startswith('data_') and 'loop_' in clean_text