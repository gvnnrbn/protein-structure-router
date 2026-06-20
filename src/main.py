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


# --- MOCK DATA GENERATOR ---
def _generate_mock_clusters() -> list:
    """
    Generates a static, determinist mock of TAPO clusters that perfectly
    matches the expected JSON contract for the frontend.
    """
    return [
        {
            "cluster_id": "cl1_selected",
            "qa_score": 0.945,
            "units": [
                {
                    "start": 10,
                    "end": 50,
                    "desc": "Cluster 1 - Unit 1 ",
                    "hex": "#ff00ff",
                    "rgb": {"r": 255, "g": 0, "b": 255}
                },
                {
                    "start": 51,
                    "end": 90,
                    "desc": "Cluster 1 - Unit 2",
                    "hex": "#cc00cc",
                    "rgb": {"r": 204, "g": 0, "b": 204}
                }
            ]
        },
        {
            "cluster_id": "cl2_selected",
            "qa_score": 0.812,
            "units": [
                {
                    "start": 150,
                    "end": 200,
                    "desc": "Cluster 2 - Unit 1",
                    "hex": "#00c8ff",
                    "rgb": {"r": 0, "g": 200, "b": 255}
                },
                {
                    "start": 201,
                    "end": 250,
                    "desc": "Cluster 2 - Unit 2",
                    "hex": "#0099cc",
                    "rgb": {"r": 0, "g": 153, "b": 204}
                }
            ]
        }
    ]

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
    # tapo_results = [_generate_mock_clusters() for _ in target_chains]
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
    # ==========================================================
    # THROWAWAY DEBUG CODE: Save the response locally to src/
    # ==========================================================
    try:
        debug_path = os.path.join("src", "debug_output.json")
        with open(debug_path, "w", encoding="utf-8") as debug_file:
            json.dump(response_payload, debug_file, indent=4)
        print(f"\n[DEBUG SUCCESS] The exact JSON response was saved locally to: {debug_path}\n")
    except Exception as e:
        print(f"[DEBUG ERROR] Could not save local debug file: {e}")
    # ==========================================================

    # Finally, return it to the browser/frontend
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