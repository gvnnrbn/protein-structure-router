# src/tapo_runner.py
import asyncio

async def run_tapo_analysis(pdb_content: str):
    """
    Receives clean PDB content and returns the array of repeats.
    """
    print("[TAPO MOCK] Processing...")
    await asyncio.sleep(1)
    
    
    return [
        { "start": 373, "end": 474, "desc": "Repeat Unit 1", "hex": "#ff00ff", "rgb": { "r": 255, "g": 0, "b": 255 } },
        { "start": 476, "end": 582, "desc": "Repeat Unit 2", "hex": "#00c8ff", "rgb": { "r": 0, "g": 200, "b": 255 } },
        { "start": 583, "end": 693, "desc": "Repeat Unit 3", "hex": "#ffa500", "rgb": { "r": 255, "g": 165, "b": 0 } }
    ]