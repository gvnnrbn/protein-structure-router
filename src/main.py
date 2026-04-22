from fastapi import FastAPI, Form, File, UploadFile, Request
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
from typing import Optional
from enum import Enum

from src.router import structure_router
from src.converters import clean_pdb_text

app = FastAPI()

class AllowedQueryTypes(str, Enum):
    """Defines the only acceptable values for query_type."""
    file = "file"
    id = "id"
    sequence = "sequence"

@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    error_messages = []
    for err in exc.errors():
        field = err["loc"][-1] if len(err["loc"]) > 0 else "unknown"
        msg = err["msg"]
        
        if field == "query_type":
            error_messages.append("The query type must be 'file', 'id', or 'sequence'.")
        elif field == "chain_id":
            error_messages.append("The chain identifier (chain_id) is required and must be a single alphanumeric character.")
        else:
            error_messages.append(f"Error in '{field}': {msg}")

    return JSONResponse(
        status_code=400,
        content={
            "status": "error",
            "message": "Invalid form data.",
            "details": error_messages
        }
    )

@app.post("/api/v1/prepare-structure")
async def process_request(
    query_type: AllowedQueryTypes = Form(...),
    chain_id: str = Form(..., min_length=1, max_length=1, pattern=r"^[A-Za-z0-9]$"),
    file_upload: Optional[UploadFile] = File(None),
    text_query: Optional[str] = Form(None)
):
    """
    Receives a structural query (File, ID, or Sequence) along with a chain ID.
    Routes it to the appropriate fetcher or converter, and prepares the PDB text for TAPO.
    """
    
    # text query cannot be empty
    if query_type in [AllowedQueryTypes.id, AllowedQueryTypes.sequence]:
        if not text_query or text_query.strip() == "":
            return JSONResponse(
                status_code=400,
                content={"status": "error", "message": "You must provide text to search by ID or Sequence."}
            )

    # file must be uploaded for file query
    if query_type == AllowedQueryTypes.file:
        if not file_upload or file_upload.filename == "":
            return JSONResponse(
                status_code=400,
                content={"status": "error", "message": "You must upload a file (.pdb, .cif, .fasta)."}
            )

    # Process the file if it exists
    file_content = None
    if file_upload and file_upload.filename != "":
        content_bytes = await file_upload.read()
        file_content = content_bytes.decode('utf-8')
    
    q_type_str = query_type.value 
    
    raw_structure = structure_router(q_type_str, text_query, file_content)
    
    if not raw_structure:
        return JSONResponse(
            status_code=404,
            content={"status": "error", "message": "Unrecognized format or structure not found in databases."}
        )
        
    # Handle Multiple Options from AlphaFold
    if isinstance(raw_structure, dict) and raw_structure.get("status") == "multiple_choices":
        return {
            "status": "multiple_choices",
            "message": "Multiple structural models found. Please choose one.",
            "chain_id": chain_id.upper(), 
            "options": raw_structure.get("options")
        }
        
    # Handle returned PDB string
    final_pdb = clean_pdb_text(raw_structure)
    
    return {
        "status": "success",
        "chain_id": chain_id.upper(),
        "pdb_preview": final_pdb[:300] + "..." 
    }