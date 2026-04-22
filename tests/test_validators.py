from src.validators import (
    is_valid_pdb_id,
    is_valid_uniprot_accession,
    is_valid_alphafold_id,
    is_fasta_format,
    is_pdb_format,
    is_mmcif_format
)

#########################################################
#       FILE FORMAT VALIDATORS
#######################################################
def test_is_fasta_format():
    valid_fasta = ">Header details\nMTEYKLVVVG"
    invalid_fasta = "MTEYKLVVVG\nNo header"
    
    assert is_fasta_format(valid_fasta) is True
    assert is_fasta_format(invalid_fasta) is False

def test_is_pdb_format():
    valid_pdb = "HEADER    STRUCTURAL PROTEIN\nATOM      1  N"
    invalid_pdb = "Just some random text\nATOM"
    
    assert is_pdb_format(valid_pdb) is True
    assert is_pdb_format(invalid_pdb) is False

def test_is_mmcif_format():
    valid_cif = "data_1XYZ\n#\nloop_\n_atom_site.group_PDB"
    
    assert is_mmcif_format(valid_cif) is True
    assert is_mmcif_format("HEADER PDB FILE") is False

#########################################################
#       ID VALIDATORS
#######################################################
def test_is_valid_pdb_id():
    # Legacy format tests
    assert is_valid_pdb_id("1LM5") is True
    assert is_valid_pdb_id("9XYZ") is True
    
    # 2028 Extended format tests
    assert is_valid_pdb_id("pdb_1000axyz") is True
    assert is_valid_pdb_id("PDB_1000AXYZ") is True  
    # Invalid format tests
    assert is_valid_pdb_id("P159") is False      
    assert is_valid_pdb_id("pdb_123") is False    
    assert is_valid_pdb_id("0ABC") is False      

def test_is_valid_uniprot_accession():
    assert is_valid_uniprot_accession("P15924") is True
    assert is_valid_uniprot_accession("p15924") is True      
    assert is_valid_uniprot_accession("A0A024RBG1") is True
    assert is_valid_uniprot_accession("P15924-2") is True   
    assert is_valid_uniprot_accession("INVALID") is False