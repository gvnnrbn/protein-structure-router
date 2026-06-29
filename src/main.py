from fastapi import FastAPI, Form, File, UploadFile, Request
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware


from src.router import structure_router
from src.tapo_runner import run_tapo_analysis
import asyncio

from pydantic import BaseModel
from typing import Optional
import requests

app = FastAPI()
origins = [
    "http://localhost:5173",
    "http://localhost:5174",
    "http://127.0.0.1:8000",    
    "http://localhost:8000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,            # Allows specific origins
    allow_credentials=True,           # Allows cookies and auth headers
    allow_methods=["*"],              # Allows all HTTP methods (GET, POST, etc.)
    allow_headers=["*"],              # Allows all request headers
)

class DetectRepeatsRequest(BaseModel):
    # Metadata
    protein_id: str
    id_type: str
    input_format: str
    
    # Specific choice properties
    choice_type: str  # "chains" or "alphafold_models"
    chain_id: str
    sequence: str
    length: int
    
    # Payload depending on choice_type
    pdb_found: Optional[str] = None
    pdb_url: Optional[str] = None

async def handle_router_result(result_dict: dict):
    """
    Shared helper function to process the router's output,
    execute TAPO if a single chain is found, and build the final response.
    """
    # If the router detected Scenario B or C (Requires user input)
    if result_dict.get("status") == "multiple_choices":
        return JSONResponse(content=result_dict)
        
    # If the router detected Scenario A (Single Chain - Auto execute)
    elif result_dict.get("status") == "success":
        pdb_content = result_dict.get("pdb_found")
        chain_id = result_dict.get("chain_id", "A")
        
        # Immediate execution of detector
        tapo_results = await run_tapo_analysis(pdb_content, chain_id, result_dict.get("protein_id", "unknown"))
        
        # Build the final payload expected by StructureView.tsx
        return JSONResponse(content={
            "status": "success",
            "input_format": result_dict.get("input_format"),
            "protein_id": result_dict.get("protein_id"),
            "id_type": result_dict.get("id_type"),
            "chain_id": chain_id,
            "sequence": result_dict.get("sequence", ""),
            "length": result_dict.get("length", 0),
            "repeats": tapo_results, 
            "pdb_found": pdb_content
        })
        
    else:
        return JSONResponse(content=result_dict, status_code=400)
    
# ---------------------------------------------------------
# ENDPOINT 1: TEXT QUERY ONLY
# ---------------------------------------------------------

@app.post("/api/prepare-structure/text")
async def process_text_request(text_query: str = Form(...)): 
    try:
        result_dict = structure_router("text", text_query=text_query)
        return await handle_router_result(result_dict)
    except Exception as e:
        print(f"[MAIN] Error processing text: {e}")
        return JSONResponse(
            content={"status": "error", "message": "Server error processing text query."}, 
            status_code=500
        )
# ---------------------------------------------------------
# ENDPOINT 2: FILE UPLOAD ONLY
# ---------------------------------------------------------

@app.post("/api/prepare-structure/file")
async def process_file_request(
    file_upload: UploadFile = File(...) 
):
    try:
        content_bytes = await file_upload.read()
        file_content = content_bytes.decode('utf-8')
        result_dict = structure_router("file", file_content=file_content)
        return await handle_router_result(result_dict)
    except Exception as e:
        print(f"[MAIN] Error processing file: {e}") 
        return JSONResponse(
            content={"status": "error", "message": "Server error processing file upload."}, 
            status_code=500
        )
# ---------------------------------------------------------
# ENDPOINT 3: DETECT REPEATS
# ---------------------------------------------------------

@app.post("/api/detect-repeats")
async def detect_repeats(request: DetectRepeatsRequest):
    try:
        pdb_content = None
        
        # SCENARIO B: Already has the raw PDB text 
        if request.choice_type == "chains":
            if not request.pdb_found:
                return JSONResponse(
                    content={"status": "error", "message": "Missing PDB content for chain selection."},
                    status_code=400
                )
            pdb_content = request.pdb_found
            
        # SCENARIO C: AlphaFold selection. The frontend only has the URL.
        elif request.choice_type == "alphafold_models":
            if not request.pdb_url:
                return JSONResponse(
                    content={"status": "error", "message": "Missing PDB URL for AlphaFold selection."},
                    status_code=400
                )
            
            print(f"[MAIN] Downloading specific AlphaFold model from: {request.pdb_url}")
            try:
                response = requests.get(request.pdb_url, timeout=15)
                if response.status_code == 200:
                    pdb_content = response.text
                else:
                    return JSONResponse(
                        content={"status": "error", "message": f"Failed to download AlphaFold model (HTTP {response.status_code})."},
                        status_code=400
                    )
            except requests.exceptions.RequestException as e:
                return JSONResponse(
                    content={"status": "error", "message": f"Network error downloading AlphaFold model: {str(e)}"},
                    status_code=500
                )
        else:
            return JSONResponse(
                content={"status": "error", "message": "Invalid choice_type provided."},
                status_code=400
            )

        print(f"[MAIN] Executing TAPO for chain {request.chain_id}...")
        
        # Execute detector
        tapo_results = await run_tapo_analysis(pdb_content, request.chain_id, request.protein_id)
        
        # Build the final payload expected by StructureView.tsx
        return JSONResponse(content={
            "status": "success",
            "input_format": request.input_format,
            "protein_id": request.protein_id,
            "id_type": request.id_type,
            "chain_id": request.chain_id,
            "sequence": request.sequence,
            "length": request.length,
            "repeats": tapo_results, 
            "pdb_found": pdb_content
        })

    except Exception as e:
        print(f"[MAIN] Error executing detect_repeats: {e}")
        return JSONResponse(
            content={"status": "error", "message": "Server error during repeat detection."}, 
            status_code=500
        )