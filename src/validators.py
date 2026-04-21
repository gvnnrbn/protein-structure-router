# Identifies type of string input

import re

def is_valid_pdb_id(pdb_id: str) -> bool:
    """
    Format: 4 alphanumeric characters. 
    """
    pattern = r'^[0-9][a-zA-Z0-9]{3}$'
    return bool(re.match(pattern, pdb_id))

def is_valid_uniprot_accession(accession: str) -> bool:
    """
    Format: 6 or 10 alphanumeric characters following specific patterns.
    """
    pattern = r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$'
    return bool(re.match(pattern, accession))

def is_valid_alphafold_id(af_id: str) -> bool:
    """
    Format: AF-<UniProtAC>-F<ModelVersion>
    """
    pattern = r'^AF-[A-Z0-9]{6,10}-F[0-9]+$'
    return bool(re.match(pattern, af_id))

def is_amino_acid_sequence(text: str) -> bool:
    """
    Checks if a string consists only of valid standard amino acid codes.
    """
    text = text.upper().strip()
    if len(text) < 20: 
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