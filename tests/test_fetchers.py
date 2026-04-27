from unittest import result

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
    result = search_rcsb_by_uniprot_id(uniprot_id)
    
    assert result["status"] == "success"
    assert result["data"] is not None
    assert "HEADER" in result["data"]
    assert "ATOM" in result["data"]

@pytest.mark.vcr()
def test_rcsb_download_by_pdb_id_valid():
    """Test PDB file download from RCSB using a valid PDB ID."""
    pdb_id = "1LM5"
    result = search_rcsb_by_pdb_id(pdb_id)
    
    assert result["status"] == "success"
    assert result["data"] is not None
    assert "HEADER" in result["data"]
    assert pdb_id in result["data"]

@pytest.mark.vcr()
def test_alphafold_multiple_choices_uniprot():
    """Test that a UniProt ID with multiple fragments returns a dictionary of choices."""
    uniprot_id = "Q5VSL9" 
    
    result = fetch_alphafold_model(uniprot_id=uniprot_id)
    
    assert result["status"] == "multiple_choices"
    assert isinstance(result["data"], list)
    assert len(result["data"]) > 1

@pytest.mark.vcr()
def test_alphafold_specific_fragment_download():
    """Test that requesting a specific AF ID downloads the exact PDB file."""
    specific_id = "AF-Q5VSL9-F1" 
    
    result = fetch_alphafold_model(specific_af_id=specific_id)
    
    assert result["status"] == "success"
    assert isinstance(result["data"], str)
    assert "ATOM" in result["data"]
    assert "toolUsed" not in result["data"]

@pytest.mark.vcr()
def test_search_by_sequence_valid_hras():
    """Test sequence search using the GTPase HRas valid sequence."""
    hras_seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPARTVETRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVIS"
    result = search_rcsb_by_sequence(hras_seq)
    
    assert result["status"] == "success"
    assert result["data"] is not None
    assert "ATOM" in result["data"]


###########################################
#       ERROR HANDLING - NEGATIVE CASES
###########################################
def test_alphafold_invalid_format_handling():
    """Test AlphaFold fetcher with invalid ID formats."""
    result = fetch_alphafold_model(uniprot_id="INVALID_ID_123")
    assert result["status"] == "error"
    assert result["data"] is None

@pytest.mark.vcr()
def test_alphafold_valid_id_no_structure():
    """Test AlphaFold fetcher with a valid-format ID that doesn't exist."""
    result = fetch_alphafold_model(specific_af_id="AF-P15924-F1")
    assert result["status"] == "error"
    assert result["data"] is None

@pytest.mark.vcr()
def test_rcsb_download_by_pdb_id_not_found():
    """Test RCSB fetcher with a non-existent PDB ID."""
    result = search_rcsb_by_pdb_id("0XXX")
    assert result["status"] == "error"
    assert result["data"] is None

@pytest.mark.vcr()
def test_rcsb_search_uniprot_no_mapping():
    """Test RCSB fetcher with a UniProt ID that has no experimental PDBs."""
    result = search_rcsb_by_uniprot_id("AX024RBG1")
    assert result["status"] == "error"
    assert result["data"] is None

@pytest.mark.vcr()
def test_search_by_sequence_no_match():
    """Test RCSB sequence fetcher with an impossible sequence."""
    weird_seq = "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"
    result = search_rcsb_by_sequence(weird_seq, identity_cutoff=0.99)
    assert result["status"] == "error"
    assert result["data"] is None

def test_sequence_validator_junk_input():
    """Ensure random junk strings fail the amino acid sequence validator."""
    assert is_amino_acid_sequence("ASDFGHJKL_NOT_A_SEQUENCE_123!") is False