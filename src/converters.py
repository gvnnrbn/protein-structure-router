import io
import re
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBParser

def parse_fasta_to_sequence(fasta_content: str) -> str:
    lines = fasta_content.strip().split('\n')
    sequence_blocks = []
    
    for line in lines:
        clean_line = line.strip()
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

def extract_pdb_metadata(result: dict, chain_id: str) -> dict:
    """
    Función auxiliar para extraer el ID, tipo de ID, secuencia y longitud de un PDB.
    Se apoya de los validadores para reconocer el formato del ID extraído.
    """
    if not result or "data" not in result or not result["data"]:
        return result
        
    pdb_content = result["data"]
    
    # Validamos que de verdad sea un formato PDB antes de iterar
    if not is_pdb_format(pdb_content):
        return result
        
    # 1. Extraer ID y su tipo (PDB o UniProt)
    protein_id = None
    id_type = "unknown"
    
    lines = pdb_content.splitlines()
    # Iteramos sobre las primeras 100 líneas del PDB donde usualmente habitan los metadatos
    for line in lines[:100]: 
        if line.startswith("HEADER"):
            parts = line.split()
            if parts:
                candidate = parts[-1]
                if is_valid_pdb_id(candidate):
                    protein_id = candidate
                    id_type = "pdbid"
                    break
        elif line.startswith("DBREF"):
            # DBREF suele tener el PDBid al inicio y la referencia Uniprot (UNP) en las columnas siguientes
            parts = line.split()
            for part in parts:
                if is_valid_uniprot_accession(part):
                    protein_id = part
                    id_type = "uniprot"
                    break
        elif line.startswith("TITLE") or line.startswith("COMPND"):
            # Respaldo por si es un archivo de AlphaFold donde agregan el Uniprot en el Title/Compnd
            for word in line.split():
                # Quitamos signos de puntuación pegados a la palabra por si acaso
                clean_word = ''.join(c for c in word if c.isalnum() or c == '-')
                if is_valid_uniprot_accession(clean_word):
                    protein_id = clean_word
                    id_type = "uniprot"
                    break
                elif is_valid_pdb_id(clean_word):
                    protein_id = clean_word
                    id_type = "pdbid"
                    break
                    
        if protein_id:
            break

    # 2. Extraer Secuencia basado en la cadena seleccionada (chain_id)
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    sequence = []
    
    # Intentamos primero con SEQRES (la forma oficial y recomendada)
    for line in lines:
        if line.startswith("SEQRES"):
            parts = line.split()
            # La columna 3 de un SEQRES contiene el identificador de cadena
            if len(parts) > 2 and parts[2] == chain_id:
                for res in parts[4:]:
                    sequence.append(three_to_one.get(res, 'X'))
                    
    # Si SEQRES no existe, recorremos los ATOMs como respaldo (Fallback)
    if not sequence:
        last_res_seq = None
        for line in lines:
            if line.startswith("ATOM"):
                # En el formato PDB, la cadena se encuentra en la columna 22 (índice 21)
                chain = line[21]
                if chain == chain_id:
                    res_name = line[17:20].strip()
                    res_seq = line[22:26].strip()
                    
                    # Evitamos agregar el mismo aminoácido repetidas veces por cada uno de sus átomos
                    if res_seq != last_res_seq:
                        sequence.append(three_to_one.get(res_name, 'X'))
                        last_res_seq = res_seq
                        
    seq_str = "".join(sequence)
    
    # 3. Asignar los nuevos valores al resultado
    result["protein_id"] = protein_id
    result["id_type"] = id_type
    result["sequence"] = seq_str
    result["length"] = len(seq_str)
    
    return result
