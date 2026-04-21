from src.converters import (
    clean_pdb_text,
)

def test_pdb_format_cleaner():
    dirty_pdb = (
        "\\nATOM      1  N   SER A2616      26.160  -7.962  45.517  1.00 36.86           N\n"
        "\\ \\nATOM      2  CA  SER A2616      27.344  -8.170  44.641  1.00 35.30           C\n"
        "\\ \\nATOM      3  C   SER A2616      27.115  -7.607  43.248  1.00 33.78           C\n"
        "\\  ADDITIONAL TEXT\n"
        "\\ \\nTER       4      SER A2616\n"
    )
    
    clean_pdb = clean_pdb_text(dirty_pdb)
    
    assert "ADDITIONAL TEXT" not in clean_pdb
    assert "\\" not in clean_pdb
    
    lines = clean_pdb.strip().split('\n')
    assert lines[0].startswith("ATOM")
    assert lines[1].startswith("ATOM")
    assert lines[-1].startswith("TER")
    
    assert len(lines[0]) > 50
