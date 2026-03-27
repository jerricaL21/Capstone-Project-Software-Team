import subprocess
import json
import glob
import os
import time
import requests
from io import StringIO
from Bio import PDB, SeqIO
from Bio.Align import PairwiseAligner

def identify_calmodulin_chain(pdb_path):

    # ── Convert PDB to FASTA using SeqIO pdb-atom parser ─────────────────────
    records = list(SeqIO.parse(pdb_path, "pdb-atom")) # Converts the PDB file directly to FASTA format using BioPython's built-in pdb-atom parser. Returns one record per chain.

    if not records: # Stops execution if no chains were found in the file.
        raise ValueError("No protein chains found in PDB file.")

    chain_seqs = {} # Empty dictionary to store chain ID 
    for record in records: # loops over each chain found in the PDB file.
        chain_id = record.id.split(":")[-1] if ":" in record.id else record.id # record.id is typically "????:A" (PDB ID + chain), extract chain letter
        seq = str(record.seq).replace("X", "")
        if seq:
            chain_seqs[chain_id] = seq  # remove unknown residues

    # ── Fetch calmodulin reference from UniProt in FASTA format (P0DP23 = human CaM-1) ───────
    url = "https://rest.uniprot.org/uniprotkb/P0DP23.fasta"
    response = requests.get(url, timeout=15)
    if not response.ok:
        raise RuntimeError("Failed to download calmodulin reference from UniProt. Check your internet connection.")
    cam_record = next(SeqIO.parse(StringIO(response.text), "fasta")) # Parses the downloaded FASTA text using SeqIO. StringIO wraps the text so SeqIO can read it like a file. 
    cam_ref = str(cam_record.seq) # Extracts the calmodulin sequence as a plain string.    

    # ── Align each chain against calmodulin reference ─────────────────────────
    aligner = PairwiseAligner()
    aligner.mode = "global"

    results = {}
    for chain_id, seq in chain_seqs.items():

        # line up the sequences first
        best = aligner.align(seq, cam_ref)[0]

        # count matches
        identical = 0
        for (qs, qe), (rs, re) in zip(best.aligned[0], best.aligned[1]):
            for q_aa, r_aa in zip(seq[qs:qe], cam_ref[rs:re]):
                if q_aa == r_aa:
                    identical += 1

        # score = how many matched out of total
        identity = (identical / len(cam_ref)) * 100

        results[chain_id] = {
            "sequence":      seq, # the full amino acid sequence of this chain
            "length":        len(seq), # how long the sequence is (number of amino acids)
            "identity":      round(identity, 1), # the match % rounded to 1 decimal e.g. 87.2
            "is_calmodulin": identity >= 80 # True if match is 80% or higher, False if not
        }

    return results


def fix_pdb_protonation(pdb_path):
    fixed_path = pdb_path.replace(".pdb", "_fixed.pdb")

    with open(pdb_path, 'r') as f:
        lines = f.readlines()

    fixed_lines = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            res_name = line[17:20].strip()
            if res_name in ["HOH", "WAT"]:
                continue
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
                chain_id = line[21].strip()
                if chain_id and chain_id not in chains:
                    chains.append(chain_id)
    return chains


def run_prppi(pdb_path, cutoff=5.0, groups=None):
    chains = get_chains(pdb_path)
    print(f"Detected chains: {chains}")

    if len(chains) < 2:
        raise ValueError(f"PDB file needs at least 2 chains, only found: {chains}")

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

    if "nano" not in stderr and process.returncode != 0:
        raise RuntimeError(f"prppi failed:\n{stderr}")

    time.sleep(2)

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


def get_residue_numbers(pdb_path):
    residues = {}
    seen     = {}

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                chain_id = line[21].strip()
                resnum   = line[22:26].strip()
                key      = f"{chain_id}{resnum}"

                if chain_id not in seen:
                    seen[chain_id]     = set()
                    residues[chain_id] = []

                if key not in seen[chain_id]:
                    seen[chain_id].add(key)
                    residues[chain_id].append(key)

    return residues


def fix_pdb_protonation(pdb_path): # This function cleans and standardises a PDB file before passing it to OSPREY. 
    fixed_path = pdb_path.replace(".pdb", "_fixed.pdb") # Creates a new output filename
    
    with open(pdb_path, 'r') as f: # Loads the entire file into memory as a list of strings, one per line.
        lines = f.readlines()
    
    fixed_lines = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"): # Only processes coordinate lines — skips REMARK, SEQRES, HEADER, CONECT, and all other record types entirely. Those pass through unchanged.
            res_name = line[17:20].strip() # Columns 17–20 in a PDB line are the residue name. 
            # Skip water only
            if res_name in ["HOH", "WAT"]: # If it's water, continue skips appending it to fixed_lines — so it simply disappears from the output file.
                continue
            # Fix histidine protonation
            if res_name in ["HIS", "HISA", "HISD", "HISE", "HISH"]: #  every histidine line get overwritten with HIE regardless of what they were before
                line = line[:17] + "HIE" + line[20:]
        fixed_lines.append(line)
    
    with open(fixed_path, 'w') as f: # Writes the cleaned file
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

def get_residue_numbers(pdb_path):
    
    residues = {}
    seen     = {}

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                chain_id = line[21].strip()
                resnum   = line[22:26].strip()  # columns 22-26 are residue sequence number
                key      = f"{chain_id}{resnum}"

                if chain_id not in seen:
                    seen[chain_id]     = set()
                    residues[chain_id] = []

                if key not in seen[chain_id]:
                    seen[chain_id].add(key)
                    residues[chain_id].append(key)

    return residues