import pytest
from unittest.mock import patch
from src.router import smart_structure_router

###########################################
#       ROUTER INTEGRATION TESTS
###########################################

# We patch the functions exactly where they are imported and used (in router.py)
@patch('src.router.search_rcsb_by_pdb_id')
def test_router_directs_pdb_id(mock_rcsb):
    """Test that a 4-character ID routes to RCSB PDB search."""
    # 1. Arrange
    mock_rcsb.return_value = "MOCKED_PDB_FILE"
    
    # 2. Act
    result = smart_structure_router("1LM5")
    
    # 3. Assert
    mock_rcsb.assert_called_once_with("1LM5") # Ensures the router passed the right ID
    assert result == "MOCKED_PDB_FILE"

@patch('src.router.search_alphafold_by_uniprot_id')
def test_router_directs_uniprot_id(mock_alphafold):
    """Test that a UniProt ID routes to AlphaFold search."""
    mock_alphafold.return_value = "MOCKED_ALPHAFOLD_FILE"
    
    result = smart_structure_router("P15924")
    
    mock_alphafold.assert_called_once_with("P15924")
    assert result == "MOCKED_ALPHAFOLD_FILE"

@patch('src.router.search_alphafold_by_uniprot_id')
def test_router_extracts_and_directs_alphafold_id(mock_alphafold):
    """Test that an AF ID is parsed correctly and routed to AlphaFold."""
    mock_alphafold.return_value = "MOCKED_ALPHAFOLD_FILE"
    
    # The user inputs the full AF ID
    result = smart_structure_router("AF-P15924-F1")
    
    # The router MUST strip the prefix/suffix and only pass "P15924"
    mock_alphafold.assert_called_once_with("P15924")
    assert result == "MOCKED_ALPHAFOLD_FILE"

@patch('src.router.search_rcsb_by_sequence')
def test_router_directs_raw_sequence(mock_sequence_search):
    """Test that a long string of valid amino acids routes to sequence search."""
    mock_sequence_search.return_value = "MOCKED_SEQUENCE_HIT"
    
    sequence = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLD"
    result = smart_structure_router(sequence)
    
    mock_sequence_search.assert_called_once_with(sequence)
    assert result == "MOCKED_SEQUENCE_HIT"

@patch('src.router.search_rcsb_by_pdb_id')
@patch('src.router.search_alphafold_by_uniprot_id')
@patch('src.router.search_rcsb_by_sequence')
def test_router_handles_invalid_input(mock_seq, mock_af, mock_pdb):
    """Test that a keyboard smash returns None and calls NO APIs."""
    
    junk_input = "!!!INVALID_JUNK_123!!!"
    result = smart_structure_router(junk_input)
    
    # 1. Assert it returns None gracefully
    assert result is None
    
    # 2. Assert that NO external API functions were triggered
    mock_pdb.assert_not_called()
    mock_af.assert_not_called()
    mock_seq.assert_not_called()