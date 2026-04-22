import pytest
from src.fetchers import (
    fetch_alphafold_model,
    search_rcsb_by_pdb_id, 
    search_rcsb_by_uniprot_id,
    search_rcsb_by_sequence,
)
from src.validators import is_amino_acid_sequence

###########################################
#       POSITIVE CASES
###########################################
@pytest.mark.vcr()
def test_rcsb_download_by_uniprot_id_valid():
    """Test PDB file download from RCSB using a valid UniProt ID."""
    uniprot_id = "P15924"
    pdb_content = search_rcsb_by_uniprot_id(uniprot_id)
    
    assert pdb_content is not None
    assert "HEADER" in pdb_content
    assert "ATOM" in pdb_content

@pytest.mark.vcr()
def test_rcsb_download_by_pdb_id_valid():
    """Test PDB file download from RCSB using a valid PDB ID."""
    pdb_id = "1LM5"
    result = search_rcsb_by_pdb_id(pdb_id)
    
    assert result is not None
    assert "HEADER" in result
    assert pdb_id in result

@pytest.mark.vcr()
def test_alphafold_multiple_choices_uniprot():
    """Test that a UniProt ID with multiple fragments returns a dictionary of choices."""
    uniprot_id = "Q5VSL9" 
    
    result = fetch_alphafold_model(uniprot_id=uniprot_id)
    
    # Must return a dictionary, NOT a PDB string
    assert isinstance(result, dict)
    assert result.get("status") == "multiple_choices"
    assert len(result.get("options")) > 1

@pytest.mark.vcr()
def test_alphafold_specific_fragment_download():
    """Test that requesting a specific AF ID downloads the exact PDB file."""
    specific_id = "AF-Q5VSL9-F1" 
    
    result = fetch_alphafold_model(specific_af_id=specific_id)
    
    # Must return the PDB string directly
    assert isinstance(result, str)
    assert "ATOM" in result
    assert "toolUsed" not in result # check that it's not JSON metadata

@pytest.mark.vcr()
def test_search_by_sequence_valid_hras():
    """Test sequence search using the GTPase HRas valid sequence."""
    hras_seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPARTVETRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVIS"
    result = search_rcsb_by_sequence(hras_seq)
    
    assert result is not None
    assert "ATOM" in result


###########################################
#       ERROR HANDLING - NEGATIVE CASES
###########################################
def test_alphafold_invalid_format_handling():
    """Test AlphaFold fetcher with invalid ID formats."""
    assert fetch_alphafold_model(uniprot_id="INVALID_ID_123") is None

@pytest.mark.vcr()
def test_alphafold_valid_id_no_structure():
    """Test AlphaFold fetcher with a valid-format ID that doesn't exist."""
    assert fetch_alphafold_model(specific_af_id="AF-P15924-F1") is None

@pytest.mark.vcr()
def test_rcsb_download_by_pdb_id_not_found():
    """Test RCSB fetcher with a non-existent PDB ID."""
    assert search_rcsb_by_pdb_id("0XXX") is None

@pytest.mark.vcr()
def test_rcsb_search_uniprot_no_mapping():
    """Test RCSB fetcher with a UniProt ID that has no experimental PDBs."""
    assert search_rcsb_by_uniprot_id("A0A024RBG1") is None

@pytest.mark.vcr()
def test_search_by_sequence_no_match():
    """Test RCSB sequence fetcher with an impossible sequence."""
    weird_seq = "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"
    assert search_rcsb_by_sequence(weird_seq, identity_cutoff=0.99) is None

def test_sequence_validator_junk_input():
    """Ensure random junk strings fail the amino acid sequence validator."""
    assert is_amino_acid_sequence("ASDFGHJKL_NOT_A_SEQUENCE_123!") is False