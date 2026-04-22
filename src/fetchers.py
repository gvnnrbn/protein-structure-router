# Structure downloader
from typing import Union, Dict

import requests

def search_rcsb_by_pdb_id(pdb_id):
    """
    Downloads a PDB file directly from the RCSB PDB repository.
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.text
    except requests.exceptions.RequestException as e:
        print(f"[FETCHER ERROR] Connection error while downloading PDB {pdb_id}: {e}")
    return None

def search_rcsb_by_uniprot_id(uniprot_id):
    """
    Searches for PDB entities using a UniProt Accession and downloads the first associated PDB file.
    """
    print(f"\n [FETCHER] Searching UniProt ID '{uniprot_id}' in RCSB Search API...")
    
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query_payload = {
      "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "operator": "exact_match",
              "value": uniprot_id, 
              "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
            }
          },
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "operator": "exact_match",
              "value": "UniProt",
              "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"
            }
          }
        ]
      },
      "return_type": "polymer_entity"
    }
    
    try:
        response_search = requests.post(search_url, json=query_payload, timeout=10)
        
        if response_search.status_code == 200:
            data = response_search.json()
            results = data.get("result_set", [])
            
            if not results:
                print(f"[FETCHER] Error: No PDB structures found for UniProt ID {uniprot_id}")
                return None
                
            # RCSB returns identifier as "1A00_1" (PDB Code + Entity ID). 
            # We split by underscore to keep only the 4-character code ("1A00").
            first_hit = results[0]["identifier"]
            clean_pdb_id = first_hit.split("_")[0] 
            
            print(f" Success: Mapped to PDB ID: {clean_pdb_id}. Starting download...")
            
            return search_rcsb_by_pdb_id(clean_pdb_id)
            
        else:
            print(f"[FETCHER] Error in RCSB Search API: Status Code {response_search.status_code}")
            return None
            
    except requests.exceptions.RequestException as e:
        print(f"[FETCHER] Connection Error: {e}")
        return None

# def search_alphafold_by_uniprot_id(uniprot_id):
#     """
#     Queries the AlphaFold API with a UniProt ID, extracts the pdbUrl for the exact sequence, and downloads the corresponding PDB file.
#     """
#     api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    
#     try:
#         response_api = requests.get(api_url, headers={'accept': 'application/json'}, timeout=10)
        
#         if response_api.status_code == 200:
#             data = response_api.json()
            
#             if not data:
#                 print(f"[FETCHER] Error: No AlphaFold model found for {uniprot_id}")
#                 return None
                
#             # Search for the exact entry (avoiding isoforms unless requested)
#             # Default to the first one but attempt to find an exact match. 
#             # TODO: RETURN TO USER TO PICK ISOFORM
#             selected_model = data[0] 
#             for entry in data:
#                 if entry.get("uniprotAccession") == uniprot_id:
#                     selected_model = entry
#                     break
                    
#             pdb_url = selected_model.get("pdbUrl")
            
#             if not pdb_url:
#                 print(f"[FETCHER] Error: No PDB URL found in AlphaFold metadata for {uniprot_id}")
#                 return None
                
#             # Download PDB file
#             response_pdb = requests.get(pdb_url, timeout=10)
            
#             if response_pdb.status_code == 200:
#                 print(f" Success: AlphaFold model for {uniprot_id} downloaded.")
#                 return response_pdb.text
                
#         return None
        
#     except requests.exceptions.RequestException as e:
#         print(f"[FETCHER] Network Error (AlphaFold): {e}")
#         return None

def fetch_alphafold_model(uniprot_id: str = None, specific_af_id: str = None) -> Union[str, Dict]:
    """
    Queries AlphaFold by UniProt ID.
    - If a specific AF ID is requested, it downloads that exact model.
    - If no specific AF ID is requested and there is only 1 result, it downloads it.
    - If no specific AF ID is requested and there are MULTIPLE results, it returns a dict of options.
    """
    if specific_af_id:
        id = specific_af_id    
    else:
        id = uniprot_id
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{id}"
    
    try:
        response_api = requests.get(api_url, headers={'accept': 'application/json'}, timeout=10)
        
        if response_api.status_code == 200:
            data = response_api.json()
            
            if not data:
                print(f"[FETCHER] Error: No AlphaFold model found for {id}")
                return None
            
            # SCENARIO 1: AF id    
            if specific_af_id:
                for entry in data:
                    if entry.get("modelEntityId") == specific_af_id:
                        return download_pdb_from_alphafold(entry.get("pdbUrl"), id)
                        
                print(f"[FETCHER] Error: Specific model {id} not found in AlphaFold data.")
                return None
            
            # SCENARIO 2: UniProt ID with MULTIPLE models
            if len(data) > 1:
                print(f"[FETCHER] Multiple models found for {id}. Prompting user selection.")
                # We return a dictionary indicating multiple choices, NOT a PDB string
                options = [
                    {
                        "af_id": entry.get("modelEntityId"), 
                        "uniprot": entry.get("uniprotAccession"),
                        "uniProtSequence": entry.get("uniprotSequence")
                    } for entry in data
                ]
                return {"status": "multiple_choices", "options": options}
            
            # SCENARIO 3: The user gave a UniProt ID, and there is exactly ONE model
            return download_pdb_from_alphafold(data[0].get("pdbUrl"), data[0].get("modelEntityId"))
                
        return None
        
    except requests.exceptions.RequestException as e:
        print(f"[FETCHER] Network Error (AlphaFold): {e}")
        return None
    
def download_pdb_from_alphafold(pdb_url: str, log_id: str) -> str:
    if not pdb_url:
        return None
    try:
        response_pdb = requests.get(pdb_url, timeout=10)
        if response_pdb.status_code == 200:
            print(f"[FETCHER] Success: AlphaFold model {log_id} downloaded.")
            return response_pdb.text
        return None
    except requests.exceptions.RequestException as e:
        print(f"[FETCHER] Network Error downloading PDB: {e}")
        return None
    
##########################################
#           SEQUENCE
##########################################
def search_rcsb_by_sequence(sequence: str, identity_cutoff: float = 0.9):
    """
    Performs a sequence search on RCSB PDB and downloads the top hit.
    """
    print(f"\n [FETCHER] Performing RCSB Sequence Search (Min Identity: {identity_cutoff*100}%)...")
    
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query_payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": identity_cutoff,
                "sequence_type": "protein",
                "value": sequence.strip().upper()
            }
        },
        "request_options": {
            "scoring_strategy": "sequence"
        },
        "return_type": "polymer_entity"
    }

    try:
        response = requests.post(search_url, json=query_payload, timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            results = data.get("result_set", [])
            
            if not results:
                print("[FETCHER] Error: No PDB matches found for this sequence.")
                return None
            
            # Extract top hit (RCSB sorts by score/identity by default)
            top_hit_id = results[0]["identifier"].split("_")[0]
            print(f" Success: Best match found: {top_hit_id} (Score: {results[0].get('score')})")
            
            # Reuse your direct download function
            return search_rcsb_by_pdb_id(top_hit_id)
            
        return None
    except requests.exceptions.RequestException as e:
        print(f"[FETCHER] Sequence Search Network Error: {e}")
        return None