import subprocess
import json
import glob
import os
import time

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

    cmd = [
        "prppi", pdb_path,
        "--side_1", str(side_1),
        "--side_2", str(side_2),
        "--cutoff", str(cutoff),
    ]

    if groups is not None:
        cmd += ["--groups", str(groups)]

    working_dir = os.path.dirname(os.path.abspath(pdb_path))

    process = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=working_dir
    )

    stdout, stderr = process.communicate(input="\n", timeout=60)

    # Ignore nano error on Windows — prppi still ran successfully
    if "nano" not in stderr and process.returncode != 0:
        raise RuntimeError(f"prppi failed:\n{stderr}")

    time.sleep(2)

    # Search inside the 'result' subfolder where prppi saves JSON
    result_dir = os.path.join(working_dir, "result")
    search_dir = result_dir if os.path.exists(result_dir) else working_dir

    print(f"Searching for JSON in: {search_dir}")

    json_files = glob.glob(os.path.join(search_dir, "*.json"))
    if json_files:
        latest_json = max(json_files, key=os.path.getctime)
        print(f"JSON found at: {latest_json}")
        return latest_json
    else:
        print(f"No JSON found in: {search_dir}")
        return None