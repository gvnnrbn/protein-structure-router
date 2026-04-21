import pytest
from src.fetchers import (
    search_rcsb_by_pdb_id, 
    search_rcsb_by_uniprot_id,
    search_alphafold_by_uniprot_id, 
    search_rcsb_by_sequence,
)
from src.validators import (
    is_amino_acid_sequence
)


###########################################
#       POSITIVE CASES
###########################################
@pytest.mark.vcr()
def test_rcsb_download_by_uniprot_id_valid():
    uniprot_id = "P15924"
    print(f"\n[TEST START] -> Fetching RCSB structure via UniProt ID: {uniprot_id}")
    
    pdb_content = search_rcsb_by_uniprot_id(uniprot_id)
    
    # 1. Validate content is returned
    assert pdb_content is not None, "Error: Result should not be None"
    
    # 2. Validate content is a genuine PDB file
    # This ensures the 'extract 4-char ID and download' flow worked
    assert "HEADER" in pdb_content, "Error: Standard PDB Header missing"
    assert "ATOM" in pdb_content, "Error: No atomic coordinates found in file"
    
    print(f"[TEST SUCCESS] -> RCSB structure for {uniprot_id} verified.")

@pytest.mark.vcr()
def test_rcsb_download_by_pdb_id_valid():
    pdb_id = "1LM5"
    print(f"\n[TEST START] -> Direct download from RCSB using PDB ID: {pdb_id}")
    
    result = search_rcsb_by_pdb_id(pdb_id)
    
    assert result is not None, "Error: Download failed"
    assert "HEADER" in result, "Error: Invalid PDB format"
    assert pdb_id in result, f"Error: Content does not match expected ID {pdb_id}"
    
    print(f"[TEST SUCCESS] -> PDB ID {pdb_id} downloaded and verified successfully.")


@pytest.mark.vcr()
def test_alphafold_download_by_uniprot_id_valid():
    uniprot_id = "P15924"
    print(f"\n[TEST START] -> Fetching full AlphaFold model for UniProt ID: {uniprot_id}")
    
    pdb_content = search_alphafold_by_uniprot_id(uniprot_id)
    
    # 1. Validate not null
    assert pdb_content is not None, "Error: AlphaFold download returned None"
    
    # 2. Validate it is NOT JSON (common error if API metadata is returned instead of the file)
    assert "toolUsed" not in pdb_content, "Error: Downloaded content is JSON metadata, not PDB text"
    
    # 3. Validate classic PDB structure
    assert "ATOM" in pdb_content, "Error: PDB file missing ATOM records"
    assert "REMARK" in pdb_content, "Error: PDB file missing REMARK metadata"
    
    print(f"[TEST SUCCESS] -> AlphaFold model downloaded and verified successfully.")

@pytest.mark.vcr()
def test_search_by_sequence_valid_hras():
    """Test using the GTPase HRas sequence provided in RCSB docs."""
    hras_seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPARTVETRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVIS"
    
    print("\n[TEST START] -> Searching for HRas sequence match...")
    result = search_rcsb_by_sequence(hras_seq)
    
    assert result is not None
    assert "ATOM" in result
    print("[TEST SUCCESS] -> Sequence search returned a valid PDB structure.")


###########################################
#       ERROR HANDLING & NEGATIVE CASES
###########################################

def test_alphafold_invalid_format_handling():
    """Test that the function rejects strings that aren't UniProt IDs immediately."""
    invalid_id = "INVALID_ID_123"
    print(f"\n[TEST START] -> Testing AlphaFold with invalid format: {invalid_id}")
    
    result = search_alphafold_by_uniprot_id(invalid_id)
    
    assert result is None, "Error: Function should return None for invalid ID formats"
    print(f"[TEST SUCCESS] -> AlphaFold correctly rejected invalid format.")

@pytest.mark.vcr()
def test_alphafold_valid_id_no_structure():
    """
    Test a valid UniProt ID that might not have an AlphaFold prediction.
    Using a made-up but valid-looking ID or a very new entry.
    """
    # Using an ID that is technically valid in format but doesn't exist in the DB
    non_existent_uniprot = "Z99999" 
    print(f"\n[TEST START] -> Testing AlphaFold with non-existent protein: {non_existent_uniprot}")
    
    result = search_alphafold_by_uniprot_id(non_existent_uniprot)
    
    assert result is None
    print(f"[TEST SUCCESS] -> AlphaFold correctly returned None for missing structure.")

@pytest.mark.vcr()
def test_rcsb_download_by_pdb_id_not_found():
    """Test RCSB direct download with a PDB ID that doesn't exist."""
    fake_pdb = "0XXX" # Valid format, but doesn't exist
    print(f"\n[TEST START] -> Testing RCSB direct download with non-existent PDB ID: {fake_pdb}")
    
    result = search_rcsb_by_pdb_id(fake_pdb)
    
    assert result is None
    print(f"[TEST SUCCESS] -> RCSB correctly returned None for missing PDB ID.")

@pytest.mark.vcr()
def test_rcsb_search_uniprot_no_mapping():
    """Test UniProt ID that exists but has NO experimental structures in RCSB."""
    uniprot_id = "A0A024RBG1" 
    print(f"\n[TEST START] -> Testing RCSB search for UniProt with no PDB entries: {uniprot_id}")
    
    result = search_rcsb_by_uniprot_id(uniprot_id)
    
    assert result is None
    print(f"[TEST SUCCESS] -> RCSB search correctly handled 'No Results' case.")

@pytest.mark.vcr()
def test_search_by_sequence_no_match():
    """Test a valid amino acid string that is highly unlikely to have a structural match."""
    weird_seq = "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"
    
    print("\n[TEST START] -> Searching for unlikely 'WWWW...' sequence...")
    result = search_rcsb_by_sequence(weird_seq, identity_cutoff=0.99)
    
    assert result is None
    print("[TEST SUCCESS] -> Correctly returned None for no match.")

def test_sequence_validator_junk_input():
    """Ensure keyboard smashes are caught before hitting the API."""
    junk = "ASDFGHJKL_NOT_A_SEQUENCE_123!"
    assert is_amino_acid_sequence(junk) is False