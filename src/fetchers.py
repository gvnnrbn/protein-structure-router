# Structure downloader
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

def search_alphafold_by_uniprot_id(uniprot_id):
    """
    Queries the AlphaFold API, extracts the pdbUrl for the exact sequence,
    and downloads the corresponding PDB file.
    """
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    
    try:
        response_api = requests.get(api_url, headers={'accept': 'application/json'}, timeout=10)
        
        if response_api.status_code == 200:
            data = response_api.json()
            
            if not data:
                print(f"[FETCHER] Error: No AlphaFold model found for {uniprot_id}")
                return None
                
            # Step 2: Search for the exact entry (avoiding isoforms unless requested)
            # Default to the first one but attempt to find an exact match. TODO: HOW TO PICK ISOFORM
            selected_model = data[0] 
            for entry in data:
                if entry.get("uniprotAccession") == uniprot_id:
                    selected_model = entry
                    break
                    
            pdb_url = selected_model.get("pdbUrl")
            
            if not pdb_url:
                print(f"[FETCHER] Error: No PDB URL found in AlphaFold metadata for {uniprot_id}")
                return None
                
            # Step 3: Download the actual PDB file from the extracted URL
            response_pdb = requests.get(pdb_url, timeout=10)
            
            if response_pdb.status_code == 200:
                print(f" Success: AlphaFold model for {uniprot_id} downloaded.")
                return response_pdb.text
                
        return None
        
    except requests.exceptions.RequestException as e:
        print(f"[FETCHER] Network Error (AlphaFold): {e}")
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
    
    # Using the exact payload structure from RCSB documentation
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