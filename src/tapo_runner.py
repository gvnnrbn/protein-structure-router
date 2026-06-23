import asyncio
import os
import shutil
import tempfile

TAPO_DIR = os.environ.get("TAPO_DIR", "/opt/tapo")
TAPO_TMP_DIR = "/home/sgeadmin/work/tmp"
TAPO_DSSP_BANK = "/home/sgeadmin/work/bank/dssp"

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
    color_palette = [
        ("#F57600", {"r": 245, "g": 118, "b": 0}),
        ("#89CE00", {"r": 137, "g": 206, "b": 0}),
        ("#FF007D", {"r": 255, "g": 0,   "b": 125}),
        ("#F5CC00", {"r": 245, "g": 204, "b": 0}),
        ("#0174E6", {"r": 1,   "g": 116, "b": 230}),
        ("#B3C7F7", {"r": 179, "g": 199, "b": 247})
    ]
    
    # cycles through the palette if repetitions exceed its size
    selected_color = color_palette[unit_index % len(color_palette)]
    
    return selected_color[0], selected_color[1]


def parse_tapo_text(raw_text: str) -> list:
    """
    Parses the raw text output from the TAPO .o file into the JSON structure
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
                    "desc": f"Unit {unit_idx + 1}",
                    "hex": hex_c,
                    "rgb": rgb_c
                })
                
            if units:
                clusters_data.extend(units)
        
    return clusters_data


def _provide_dssp(project: str, chain: str) -> None:
    """
    Copy a pre-computed DSSP file from the bank into the TAPO TMP_DIR so that
    TAPO can find it at {TMP_DIR}/{project}{chain}.dssp.
    DSSP bank layout: {DSSP_BANK}/{pdbid[1:3]}/{pdbid}{chain}.dssp
    """
    dssp_target = os.path.join(TAPO_TMP_DIR, f"{project}{chain}.dssp")
    if os.path.exists(dssp_target):
        return

    if len(project) < 3:
        return

    mid2 = project[1:3]
    bank_path = os.path.join(TAPO_DSSP_BANK, mid2, f"{project}{chain}.dssp")
    if os.path.exists(bank_path):
        os.makedirs(TAPO_TMP_DIR, exist_ok=True)
        shutil.copy2(bank_path, dssp_target)
        print(f"[DSSP] Copied from bank: {bank_path} → {dssp_target}")


async def run_tapo_analysis(pdb_content: str, target_chain: str, protein_id: str, timeout_seconds: int = 240):
    print(f"[TAPO RUNNER] Preparing direct Java execution for chain {target_chain}...")

    # TAPO uses project name to look up DSSP at {TMP_DIR}/{project}{chain}.dssp.
    # The bank stores files as {pdbid_lower}{chain}.dssp, so we use lowercase PDB ID.
    project = protein_id.lower()

    with tempfile.TemporaryDirectory() as temp_dir:
        input_pdb_path = os.path.join(temp_dir, f"{project}_{target_chain}.pdb")
        output_path = os.path.join(temp_dir, f"{project}_{target_chain}.o")
        with open(input_pdb_path, "w", encoding="utf-8") as f:
            f.write(pdb_content)
        # Provide DSSP from bank when available (dramatically improves scoring accuracy)
        _provide_dssp(project, target_chain)

        classpath = "target/*:dependencies/bugs/vecmath-1.3.1.jar:dependencies/*"
        
        output_file_name = f"{protein_id}_{target_chain}.o"

        cmd = [
            "java",
            "-Xmx2042m",
            "-Djava.awt.headless=true",
            "-Dlog4j.skipJansi=true",
            "-cp", classpath,
            "apps.TaPo",
            "-f", input_pdb_path,
            "-p", project,
            "-c", target_chain,
            "-nCore", "1",
            "-o", output_path,
        ]

        print(f"[TAPO RUNNER] Command: java ... apps.TaPo -f {input_pdb_path} -p {project} -c {target_chain} -o {output_path}")

        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=TAPO_DIR,
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