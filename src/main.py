from fastapi import FastAPI, Form, File, UploadFile, Request
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware


from src.router import structure_router
from src.converters import clean_pdb_text, get_chains_from_pdb
from src.tapo_runner import run_tapo_analysis
import asyncio

import json
import os

app = FastAPI()

_frontend_url = os.environ.get("FRONTEND_URL")
origins = [_frontend_url] if _frontend_url else ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
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
    
    # 1. Determine which physical chains to process
    target_chains = []
    if chain_id.upper() == "ALL":
        target_chains = get_chains_from_pdb(final_pdb)
        # Fallback in case BioPython fails to find chains
        if not target_chains:
            target_chains = ["A"] 
    else:
        target_chains = [chain_id.upper()]
        
    # 2. ASYNC ORCHESTRATION: Trigger TAPO containers for all target chains simultaneously
    tasks = [run_tapo_analysis(final_pdb, ch,result["protein_id"]) for ch in target_chains]
    tapo_results = await asyncio.gather(*tasks)

    # 3. Aggregation: Map the output JSON from each container to its respective chain ID
    repeats_by_chain = {ch: res for ch, res in zip(target_chains, tapo_results)}
    
    # 4. Construct Final Response Payload
    if chain_id.upper() != "ALL":
        response_payload = {
            "status": "success",
            "input_format": result["format"],
            "protein_id": result["protein_id"],
            "id_type": result["id_type"],
            "chain_id": chain_id.upper(),
            "sequence": result["sequence"],
            "length": result["length"],
            "repeats": repeats_by_chain[chain_id.upper()],
            "pdb_found": final_pdb
        }
    else:
        response_payload = {
            "status": "success",
            "input_format": result["format"],
            "protein_id": result["protein_id"],
            "id_type": result["id_type"],
            "length": result["length"],
            "repeats": repeats_by_chain[chain_id.upper()],
            "chain_id": "ALL",
            "chains_data": repeats_by_chain, 
            "pdb_found": final_pdb
        }
    return response_payload
        

# ---------------------------------------------------------
# ENDPOINT 1: TEXT QUERY ONLY
# ---------------------------------------------------------
@app.post("/api/prepare-structure/text")
async def process_text_request(
    chain_id: str = Form(..., pattern=r"^([A-Za-z0-9]|[Aa][Ll][Ll])$"),
    text_query: str = Form(...) 
):
    """Handles ID and Sequence queries."""
    result_dict = structure_router("text", text_query=text_query, chain_id=chain_id)
    
    return await _handle_router_response(result_dict, chain_id)

# ---------------------------------------------------------
# ENDPOINT 2: FILE UPLOAD ONLY
# ---------------------------------------------------------
@app.post("/api/prepare-structure/file")
async def process_file_request(
    chain_id: str = Form(..., pattern=r"^([A-Za-z0-9]|[Aa][Ll][Ll])$"),
    file_upload: UploadFile = File(...) 
):
    """Handles PDB, mmCIF, and FASTA file uploads."""
    content_bytes = await file_upload.read()
    file_content = content_bytes.decode('utf-8')
    
    result_dict = structure_router("file", file_content=file_content, chain_id=chain_id)
    
    return await _handle_router_response(result_dict, chain_id)