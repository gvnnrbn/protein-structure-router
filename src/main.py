from src.router import smart_structure_router
from src.utils import clean_pdb_text

def process_request(user_input: str):
    
    raw_structure = smart_structure_router(user_input)
    
    if not raw_structure:
        return "Error: Formato no reconocido o estructura no encontrada."
        
    final_pdb = clean_pdb_text(raw_structure)
    
    # return send_to_detector(final_pdb)