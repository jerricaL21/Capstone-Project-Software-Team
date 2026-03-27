import subprocess
import json
import glob
import os
import time
import requests
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner

def identify_calmodulin_chain(pdb_path):
    """
    Runs on the raw uploaded PDB file (before protonation fix).
    Extracts each chain's amino acid sequence using Biopython PPBuilder,
    then aligns each against the human calmodulin reference (UniProt P0DP23).
    Returns a dict keyed by chain ID with identity scores and classification.
    """
    # ── Extract sequences from each chain ────────────────────────────────────
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_path)
    builder = PPBuilder()

    chain_seqs = {}
    for model in structure:
        for chain in model:
            # PPBuilder automatically ignores waters and HETATM records
            polypeptides = builder.build_peptides(chain)
            if polypeptides:
                full_seq = "".join(str(pp.get_sequence()) for pp in polypeptides)
                chain_seqs[chain.id] = full_seq

    if not chain_seqs:
        raise ValueError("No protein chains found in PDB file.")

    # ── Fetch calmodulin reference from UniProt (P0DP23 = human CaM-1) ───────
    url = "https://rest.uniprot.org/uniprotkb/P0DP23.fasta"
    response = requests.get(url, timeout=15)
    if not response.ok:
        raise RuntimeError("Failed to download calmodulin reference from UniProt. Check your internet connection.")
    fasta_lines = response.text.strip().split("\n")
    cam_ref = "".join(fasta_lines[1:])  # skip the >header line

    # ── Configure aligner (equivalent to globalms 2, -1, -10, -0.5) ─────────
    aligner = PairwiseAligner()
    aligner.mode            = "global"
    aligner.match_score     = 2
    aligner.mismatch_score  = -1
    aligner.open_gap_score  = -10
    aligner.extend_gap_score = -0.5

    # ── Align each chain against calmodulin reference ─────────────────────────
    results = {}
    for chain_id, seq in chain_seqs.items():
        alignments = aligner.align(seq, cam_ref)

        if alignments:
            best = alignments[0]
            aligned_query, aligned_ref = best.aligned

            # Count identical positions
            identical = 0
            align_len = 0
            for (qs, qe), (rs, re) in zip(aligned_query, aligned_ref):
                for q_aa, r_aa in zip(seq[qs:qe], cam_ref[rs:re]):
                    align_len += 1
                    if q_aa == r_aa:
                        identical += 1

            identity = (identical / align_len) * 100 if align_len > 0 else 0.0
        else:
            identity = 0.0

        results[chain_id] = {
            "sequence":      seq,
            "length":        len(seq),
            "identity":      round(identity, 1),
            "is_calmodulin": identity >= 80
        }

    return results

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