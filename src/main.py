from fastapi import FastAPI, Form, File, UploadFile, Request
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware


from src.router import structure_router
from src.converters import clean_pdb_text
from src.tapo_runner import run_tapo_analysis

app = FastAPI()
origins = [
    "http://localhost:5173",
    "http://localhost:5174",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,            # Allows specific origins
    allow_credentials=True,           # Allows cookies and auth headers
    allow_methods=["*"],              # Allows all HTTP methods (GET, POST, etc.)
    allow_headers=["*"],              # Allows all request headers
)

@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    error_messages = []
    for err in exc.errors():
        field = err["loc"][-1] if len(err["loc"]) > 0 else "unknown"
        msg = err["msg"]
        
        if field == "chain_id":
            error_messages.append("The chain identifier (chain_id) is required and must be a single alphanumeric character.")
        elif field == "file_upload":
            error_messages.append("You must upload a valid file.")
        elif field == "text_query":
            error_messages.append("You must provide text to search.")
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

async def _handle_router_response(result: dict, chain_id: str):
    """Helper function to DRY up the response logic for both endpoints."""
    if result["status"] == "error":
        return JSONResponse(
            status_code=404,
            content={
                "status": "error",
                "input_format": result["format"],
                "message": result["message"]
            }
        )
        
    # Handle Multiple Options from AlphaFold
    if result["status"] == "multiple_choices":
        print("[BACKEND] Multiple choices found. Returning options to frontend.")
        return {
            "status": "multiple_choices",
            "input_format": result["format"],
            "message": "Multiple models found.",
            "chain_id": chain_id.upper(), 
            "options": result["data"]
        }
        
    # Handle returned PDB string
    print("[BACKEND] Unique PDB found. Sending to TAPO...")
    final_pdb = clean_pdb_text(result["data"])

    # TAPO mock
    repeats_data = await run_tapo_analysis(final_pdb)
    print("[BACKEND] Returning detector data...")
    return {
        "status": "success",
        "input_format": result["format"],
        "chain_id": chain_id.upper(),
        "repeats": repeats_data,
        "pdb_found": final_pdb,
    }

# ---------------------------------------------------------
# ENDPOINT 1: TEXT QUERY ONLY
# ---------------------------------------------------------
@app.post("/api/prepare-structure/text")
async def process_text_request(
    chain_id: str = Form(..., min_length=1, max_length=1, pattern=r"^[A-Za-z0-9]$"),
    text_query: str = Form(...) 
):
    """Handles ID and Sequence queries."""
    result_dict = structure_router("text", text_query=text_query)
    
    return await _handle_router_response(result_dict, chain_id)

# ---------------------------------------------------------
# ENDPOINT 2: FILE UPLOAD ONLY
# ---------------------------------------------------------
@app.post("/api/prepare-structure/file")
async def process_file_request(
    chain_id: str = Form(..., min_length=1, max_length=1, pattern=r"^[A-Za-z0-9]$"),
    file_upload: UploadFile = File(...) 
):
    """Handles PDB, mmCIF, and FASTA file uploads."""
    content_bytes = await file_upload.read()
    file_content = content_bytes.decode('utf-8')
    
    result_dict = structure_router("file", file_content=file_content)
    
    return await _handle_router_response(result_dict, chain_id)