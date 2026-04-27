from unittest import result

import pytest
from unittest.mock import patch
from src.router import structure_router

#########################################################
#       TEXT QUERY TESTS (IDs & Sequences)
#######################################################
@patch('src.router.search_rcsb_by_pdb_id')
def test_router_directs_pdb_id(mock_rcsb):
    """Test that a 4-character ID routes to RCSB PDB search."""
    mock_rcsb.return_value = {"status": "success", "data": "MOCKED_PDB_FILE", "message": None}    
    result = structure_router(query_type="text", text_query="1LM5")
    
    mock_rcsb.assert_called_once_with("1LM5")
    assert result["status"] == "success"
    assert result["data"] == "MOCKED_PDB_FILE"
    assert result["format"] == "PDB ID"

@patch('src.router.fetch_alphafold_model')
def test_router_directs_uniprot_id(mock_alphafold):
    """Test that a UniProt ID routes to AlphaFold search without specific fragment."""
    mock_options = [
        {"af_id": "AF-Q5VSL9-F1", "uniprot": "Q5VSL9", "uniProtSequence": "MEPAV..."},
        {"af_id": "AF-Q5VSL9-4-F1", "uniprot": "Q5VSL9-4", "uniProtSequence": "MEPAV..."},
        {"af_id": "AF-Q5VSL9-3-F1", "uniprot": "Q5VSL9-3", "uniProtSequence": "MEPAV..."},
        {"af_id": "AF-Q5VSL9-2-F1", "uniprot": "Q5VSL9-2", "uniProtSequence": "MEPAV..."},
    ]
    
    mock_alphafold.return_value = {
        "status": "multiple_choices", 
        "data": mock_options, 
        "message": None
    }
    result = structure_router(query_type="text", text_query="Q5VSL9")    
    mock_alphafold.assert_called_once_with(uniprot_id="Q5VSL9")
    assert result["status"] == "multiple_choices"
    assert result["format"] == "Uniprot Accession"
    assert len(result["data"]) == 4
    assert result["data"][0]["af_id"] == "AF-Q5VSL9-F1"

@patch('src.router.fetch_alphafold_model')
def test_router_extracts_and_directs_alphafold_id(mock_alphafold):
    """Test that an AF ID passes both UniProt ID and specific AF ID to fetcher."""
    mock_alphafold.return_value = {"status": "success", "data": "MOCKED_ALPHAFOLD_FILE", "message": None}
    
    result = structure_router(query_type="text", text_query="AF-Q5VSL9-F1")
    
    mock_alphafold.assert_called_once_with(specific_af_id="AF-Q5VSL9-F1")
    assert result["data"] == "MOCKED_ALPHAFOLD_FILE"
    assert result["format"] == "AlphaFold ID"

@patch('src.router.search_rcsb_by_sequence')
def test_router_directs_raw_sequence(mock_sequence_search):
    """Test that a valid amino acid sequence routes to sequence search."""
    mock_sequence_search.return_value = {"status": "success", "data": "MOCKED_SEQUENCE_HIT", "message": None}
    sequence = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLD"
    
    result = structure_router(query_type="text", text_query=sequence)
    
    mock_sequence_search.assert_called_once_with(sequence)
    assert result["data"] == "MOCKED_SEQUENCE_HIT"
    assert result["format"] == "Sequence"

@patch('src.router.search_rcsb_by_pdb_id')
@patch('src.router.fetch_alphafold_model')
@patch('src.router.search_rcsb_by_sequence')
def test_router_handles_invalid_input(mock_seq, mock_af, mock_pdb):
    """Test that invalid input returns None and triggers no APIs."""
    junk_input = "!#INVALID_INPUT123!"
    
    result = structure_router(query_type="text", text_query=junk_input)
    
    assert result["status"] == "error"
    assert result["data"] is None
    assert "Unrecognized" in result["message"]
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
    assert result["status"] == "success"
    assert result["data"] == mock_pdb_content
    assert result["format"] == "PDB"

@patch('src.router.convert_mmcif_to_pdb')
def test_router_converts_mmcif_file(mock_convert):
    """Test that an mmCIF file triggers the conversion function."""
    mock_mmcif_content = "data_1XYZ\nloop_\n_atom_site.group_PDB"
    mock_convert.return_value = {
        "status": "success", 
        "data": "CONVERTED_PDB_CONTENT", 
        "message": None
    }
    
    result = structure_router(query_type="file", file_content=mock_mmcif_content)
    
    mock_convert.assert_called_once_with(mock_mmcif_content)
    assert result["status"] == "success"
    assert result["data"] == "CONVERTED_PDB_CONTENT"
    assert result["format"] == "mmCIF"

@patch('src.router.search_rcsb_by_sequence')
def test_router_extracts_fasta_and_searches(mock_sequence_search):
    """Test that a FASTA file is parsed and sent to sequence search."""
    mock_fasta_content = ">sp|P15924|DESP_HUMAN Desmoplakin\nMTEYKLVVVGAGGVGKSALTIQLIQNHF"
    mock_sequence_search.return_value = {"status": "success", "data": "PDB_FROM_FASTA_SEQUENCE", "message": None}
    
    result = structure_router(query_type="file", file_content=mock_fasta_content)
    
    mock_sequence_search.assert_called_once_with("MTEYKLVVVGAGGVGKSALTIQLIQNHF")
    assert result["status"] == "success"
    assert result["data"] == "PDB_FROM_FASTA_SEQUENCE"
    assert result["format"] == "FASTA"