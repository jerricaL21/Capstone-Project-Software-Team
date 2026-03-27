import json
import glob
import os
import sys

# Add prppi package to path so we can import it directly
PRPPI_PATH = r"C:\prppi-main"
if PRPPI_PATH not in sys.path:
    sys.path.insert(0, PRPPI_PATH)

def fix_pdb_protonation(pdb_path):
    fixed_path = pdb_path.replace(".pdb", "_fixed.pdb")
    
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            res_name = line[17:20].strip()
            # Skip water only
            if res_name in ["HOH", "WAT"]:
                continue
            # Fix histidine protonation
            if res_name in ["HIS", "HISA", "HISD", "HISE", "HISH"]:
                line = line[:17] + "HIE" + line[20:]
        fixed_lines.append(line)
    
    with open(fixed_path, 'w') as f:
        f.writelines(fixed_lines)
    
    return fixed_path

def get_chains(pdb_path):
    """Read the PDB file and return a list of unique chain IDs."""
    chains = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21].strip()  # chain ID is always at column 21
                if chain_id and chain_id not in chains:
                    chains.append(chain_id)
    return chains


def run_prppi(pdb_path, cutoff=5.0, groups=None):
    # Auto-detect chains from the PDB file
    chains = get_chains(pdb_path)
    print(f"Detected chains: {chains}")

    if len(chains) < 2:
        raise ValueError(f"PDB file needs at least 2 chains, only found: {chains}")

    # Use first two chains automatically
    side_1 = chains[0]
    side_2 = chains[1]
    print(f"Using side_1={side_1}, side_2={side_2}")

    # Import and call prppi directly — no subprocess needed
    from prppi import sort_defined_interface

    # prppi saves result to ./result/out.json relative to cwd,
    # so we change cwd to the PDB file's directory
    working_dir = os.path.dirname(os.path.abspath(pdb_path))
    original_dir = os.getcwd()

    try:
        os.chdir(working_dir)
        json_path = sort_defined_interface.run(pdb_path, side_1, side_2, cutoff, groups)
    finally:
        os.chdir(original_dir)

    if json_path and os.path.exists(json_path):
        print(f"JSON found at: {json_path}")
        return json_path

    # Fallback: search for the JSON in case path returned is relative
    result_dir = os.path.join(working_dir, "result")
    search_dir = result_dir if os.path.exists(result_dir) else working_dir
    json_files = glob.glob(os.path.join(search_dir, "*.json"))
    if json_files:
        latest_json = max(json_files, key=os.path.getctime)
        print(f"JSON found at: {latest_json}")
        return latest_json

    print(f"No JSON found in: {search_dir}")
    return None