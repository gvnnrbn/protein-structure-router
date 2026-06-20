import asyncio
import os
import tempfile

# Pre-defined palette for up to 5 distinct tandem repeat domains (clusters)
CLUSTER_PALETTES = [
    (230, 159, 0),  
    (86, 180, 233),  
    (0, 158, 115),   
    (204, 121, 167), 
    (0, 114, 178)    
]

def get_contrasting_color(base_rgb: tuple, unit_index: int):
    """
    Generates a derived color by darkening the base color based on the unit index.
    unit_index 0 = base color. unit_index 1 = 20% darker, etc.
    """
    factor = 1.0 - (unit_index * 0.20)
    # Prevents the color from becoming completely black
    factor = max(0.4, factor) 
    
    r = int(base_rgb[0] * factor)
    g = int(base_rgb[1] * factor)
    b = int(base_rgb[2] * factor)
    
    hex_color = f"#{r:02x}{g:02x}{b:02x}"
    return hex_color, {"r": r, "g": g, "b": b}


def parse_tapo_text(raw_text: str) -> list:
    """
    Parses the raw text output from the TAPO .o file into the JSON structure
    expected by the frontend.
    
    Filters rows where the 3rd column ends with '_selected'.
    """
    clusters_data = []
    cluster_count = 0
    
    # Split the raw text into individual lines
    lines = raw_text.strip().split('\n')
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith(">"):
            continue

        parts = line.split('\t')
        if len(parts) < 8:
            continue
            
        cluster_id = parts[2].strip()
        
        # Only process the best prediction for each domain
        if cluster_id.endswith("_selected"):
            base_color = CLUSTER_PALETTES[cluster_count % len(CLUSTER_PALETTES)]
            cluster_count += 1
            
            # Extract the Repeat Units regions (Column 7, index 6)
            regions_str = parts[6].strip()
            region_parts = regions_str.split(";")
            
            units = []
            for unit_idx, reg in enumerate(region_parts):
                if not reg or "-" not in reg:
                    continue
                
                start_str, end_str = reg.split("-")
                
                try:
                    start_pos = int(start_str)
                    end_pos = int(end_str)
                except ValueError:
                    continue # Skip malformed regions 
                
                # Generate a unique shade for this repeat unit
                hex_c, rgb_c = get_contrasting_color(base_color, unit_idx)
                
                units.append({
                    "start": start_pos,
                    "end": end_pos,
                    "desc": f"Cluster {cluster_count} - Unit {unit_idx + 1}",
                    "hex": hex_c,
                    "rgb": rgb_c
                })
                
            if units:
                clusters_data.extend(units)
        
    return clusters_data

async def run_tapo_analysis(pdb_content: str, target_chain: str, protein_id: str, timeout_seconds: int = 120):
    """
    Executes the TAPO Docker container asynchronously with safety limits.
    Includes Docker resource constraints and asyncio timeouts to prevent host freezes.
    """
    print(f"[TAPO RUNNER] Preparing secure Docker execution for chain {target_chain}...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        input_pdb_name = f"{protein_id}_{target_chain}.pdb"
        input_pdb_path = os.path.join(temp_dir, input_pdb_name)
        
        with open(input_pdb_path, "w", encoding="utf-8") as f:
            f.write(pdb_content)
            
        output_file_name = f"{protein_id}_{target_chain}.o"
        
        # SAFETY MEASURE 1: Docker limits
        docker_cmd = (
            f"docker run --rm "
            f"--cpus=\"1.0\" -m=\"1g\" "
            f"-v {temp_dir}:/pdbdata "
            f"-v {temp_dir}:/workdata "
            f"-w /home/sgeadmin/save/BioApps/tapo-v1.1.3 "
            f"pdoviet/tapo:v1.1.3-alpha.2 "
            f"sh bin/run apps.TaPo '-f /pdbdata/{input_pdb_name} -p TEMP -c {target_chain} -nCore 1 -o /workdata/{output_file_name}'"
        )
        
        print(f"[TAPO RUNNER] Running command: {docker_cmd}")
        
        process = await asyncio.create_subprocess_shell(
            docker_cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        
        try:
            # SAFETY MEASURE 2: Asyncio Timeout
            stdout, stderr = await asyncio.wait_for(process.communicate(), timeout=timeout_seconds)
            
        except asyncio.TimeoutError:
            print(f"[TAPO ERROR] Process for chain {target_chain} exceeded {timeout_seconds}s timeout. Terminating...")
            process.kill()
            await process.communicate()
            return []
            
        if process.returncode != 0:
            err_msg = stderr.decode().strip()
            print(f"[TAPO ERROR] Docker failed for chain {target_chain}. Details: {err_msg}")
            return []
            
        output_path = os.path.join(temp_dir, output_file_name)
        
        if not os.path.exists(output_path):
            print(f"[TAPO ERROR] Output file was not generated for chain {target_chain}")
            return []
            
        with open(output_path, "r", encoding="utf-8") as out_f:
            tapo_raw_text = out_f.read()
            
        # Call the parser to transform the raw text into the structured JSON
        parsed_clusters = parse_tapo_text(tapo_raw_text)
        
        print(f"[TAPO RUNNER] Execution finished successfully for chain {target_chain}. Found {len(parsed_clusters)} valid clusters.")
        return parsed_clusters