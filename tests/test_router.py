import pytest
from unittest.mock import patch
from src.router import structure_router

#########################################################
#       TEXT QUERY TESTS (IDs & Sequences)
#######################################################
@patch('src.router.search_rcsb_by_pdb_id')
def test_router_directs_pdb_id(mock_rcsb):
    """Test that a 4-character ID routes to RCSB PDB search."""
    mock_rcsb.return_value = "MOCKED_PDB_FILE"
    
    result = structure_router(query_type="id", text_query="1LM5")
    
    mock_rcsb.assert_called_once_with("1LM5")
    assert result == "MOCKED_PDB_FILE"

@patch('src.router.fetch_alphafold_model')
def test_router_directs_uniprot_id(mock_alphafold):
    """Test that a UniProt ID routes to AlphaFold search without specific fragment."""
    mock_alphafold.return_value = {"status": "multiple_choices"} 
    
    result = structure_router(query_type="id", text_query="Q5VSL9")
    
    mock_alphafold.assert_called_once_with(uniprot_id="Q5VSL9")
    assert result == {"status": "multiple_choices"}

@patch('src.router.fetch_alphafold_model')
def test_router_extracts_and_directs_alphafold_id(mock_alphafold):
    """Test that an AF ID passes both UniProt ID and specific AF ID to fetcher."""
    mock_alphafold.return_value = "MOCKED_ALPHAFOLD_FILE"
    
    result = structure_router(query_type="id", text_query="AF-Q5VSL9-F1")
    
    # Must call with both arguments
    mock_alphafold.assert_called_once_with(specific_af_id="AF-Q5VSL9-F1")
    assert result == "MOCKED_ALPHAFOLD_FILE"

@patch('src.router.search_rcsb_by_sequence')
def test_router_directs_raw_sequence(mock_sequence_search):
    """Test that a valid amino acid sequence routes to sequence search."""
    mock_sequence_search.return_value = "MOCKED_SEQUENCE_HIT"
    sequence = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLD"
    
    result = structure_router(query_type="sequence", text_query=sequence)
    
    mock_sequence_search.assert_called_once_with(sequence)
    assert result == "MOCKED_SEQUENCE_HIT"

@patch('src.router.search_rcsb_by_pdb_id')
@patch('src.router.fetch_alphafold_model')
@patch('src.router.search_rcsb_by_sequence')
def test_router_handles_invalid_input(mock_seq, mock_af, mock_pdb):
    """Test that invalid input returns None and triggers no APIs."""
    junk_input = "!#INVALID_INPUT123!"
    
    result = structure_router(query_type="id", text_query=junk_input)
    
    assert result is None
    mock_pdb.assert_not_called()
    mock_af.assert_not_called()
    mock_seq.assert_not_called()


###########################################################
#       FILE UPLOAD TESTS (PDB, mmCIF, FASTA)
#########################################################
def test_router_passes_through_pdb_file():
    """Test that an uploaded PDB file is returned directly."""
    mock_pdb_content = "HEADER    STRUCTURAL PROTEIN\nATOM      1  N   ALA A   1"
    
    result = structure_router(query_type="file", file_content=mock_pdb_content)
    assert result == mock_pdb_content

@patch('src.router.convert_mmcif_to_pdb')
def test_router_converts_mmcif_file(mock_convert):
    """Test that an mmCIF file triggers the conversion function."""
    mock_mmcif_content = "data_1XYZ\nloop_\n_atom_site.group_PDB"
    mock_convert.return_value = "CONVERTED_PDB_CONTENT"
    
    result = structure_router(query_type="file", file_content=mock_mmcif_content)
    
    mock_convert.assert_called_once_with(mock_mmcif_content)
    assert result == "CONVERTED_PDB_CONTENT"

@patch('src.router.search_rcsb_by_sequence')
def test_router_extracts_fasta_and_searches(mock_sequence_search):
    """Test that a FASTA file is parsed and sent to sequence search."""
    mock_fasta_content = ">sp|P15924|DESP_HUMAN Desmoplakin\nMTEYKLVVVGAGGVGKSALTIQLIQNHF"
    mock_sequence_search.return_value = "PDB_FROM_FASTA_SEQUENCE"
    
    result = structure_router(query_type="file", file_content=mock_fasta_content)
    
    mock_sequence_search.assert_called_once_with("MTEYKLVVVGAGGVGKSALTIQLIQNHF")
    assert result == "PDB_FROM_FASTA_SEQUENCE"