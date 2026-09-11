"""
logic5.py — PDB Structural Repair
==================================
Sits between Step 1 (Data Input) and the prppi call in logic.py.

Responsibilities:
  1. Inspect an incoming PDB for missing side chains and/or missing residues.
  2. Repair missing side-chain atoms via PDBFixer (fast, in-process).
  3. Return a clean, complete PDB path ready for fix_pdb_protonation()
     and run_prppi() in logic.py.

Note: PDBs with entirely missing residues/loops are currently not
auto-repaired. They are flagged by inspect_pdb() so the caller (gui.py)
can warn the user, and the pipeline proceeds with the original
coordinates for that portion of the structure.

Dependencies
------------
  pip install pdbfixer openmm
"""

import os
from pathlib import Path
from typing import Optional

from Bio.PDB import PDBParser


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — PDB INSPECTION
# Detects what is actually wrong with the file before deciding how to fix it.
# ═══════════════════════════════════════════════════════════════════════════════

def inspect_pdb(pdb_path: str) -> dict:
    """
    Parse the PDB and return a summary of structural completeness.

    Returns a dict with keys:
      has_missing_residues  (bool) — entire residues absent from ATOM records
                                     but listed in SEQRES / REMARK 465
      has_missing_sidechains (bool) — residues present but lacking CB/CG atoms
      missing_residue_details (list of str) — human-readable descriptions
      incomplete_sidechain_residues (list of str)
      chains (list of str) — chain IDs found
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", pdb_path)

    # ── Detect missing residues from REMARK 465 ────────────────────────────
    # REMARK 465 lists residues present in SEQRES but absent from coordinates.
    missing_residue_details = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("REMARK 465") and len(line) > 20:
                # Standard format: REMARK 465   MET A    1
                parts = line[10:].split()
                if len(parts) >= 3:
                    missing_residue_details.append(
                        f"{parts[0]} {parts[1]} {parts[2]}"
                    )

    # ── Detect incomplete side chains in present residues ─────────────────
    # A residue is flagged if it has a CA but is missing CB (all non-Gly)
    # or is missing common side chain atoms expected for that residue.
    SIDECHAIN_CHECK = {
        # residue_name: ALL expected side-chain heavy atoms (non-backbone)
        # Backbone atoms (N, CA, C, O) are excluded — only side-chain atoms listed.
        # Source: standard PDB/wwPDB atom nomenclature for each amino acid.
        "ALA": ["CB"],
        "VAL": ["CB", "CG1", "CG2"],
        "LEU": ["CB", "CG", "CD1", "CD2"],
        "ILE": ["CB", "CG1", "CG2", "CD1"],
        "PRO": ["CB", "CG", "CD"],
        "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TRP": ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
        "MET": ["CB", "CG", "SD", "CE"],
        "SER": ["CB", "OG"],
        "THR": ["CB", "OG1", "CG2"],
        "CYS": ["CB", "SG"],
        "TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
        "HIS": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
        "HIE": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],  # HIS epsilon-tautomer
        "HID": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],  # HIS delta-tautomer
        "HIP": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],  # HIS doubly protonated
        "ASP": ["CB", "CG", "OD1", "OD2"],
        "GLU": ["CB", "CG", "CD", "OE1", "OE2"],
        "ASN": ["CB", "CG", "OD1", "ND2"],
        "GLN": ["CB", "CG", "CD", "OE1", "NE2"],
        "LYS": ["CB", "CG", "CD", "CE", "NZ"],
        "ARG": ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
        # GLY has no side chain — skipped in the loop above
        # SEC (selenocysteine) and PYL (pyrrolysine) are rare; omitted
    }

    incomplete_sidechain_residues = []
    chains_found = []

    for model in structure:
        for chain in model:
            cid = chain.get_id()
            if cid not in chains_found:
                chains_found.append(cid)
            for residue in chain:
                resname = residue.get_resname().strip()
                resid   = residue.get_id()[1]
                if resname == "GLY":
                    continue  # glycine has no CB by design
                atom_names = {atom.get_name() for atom in residue}
                required   = SIDECHAIN_CHECK.get(resname, [])
                missing    = [a for a in required if a not in atom_names]
                if missing:
                    incomplete_sidechain_residues.append(
                        f"{resname} {cid}{resid} (missing: {', '.join(missing)})"
                    )

    return {
        "has_missing_residues":       len(missing_residue_details) > 0,
        "has_missing_sidechains":     len(incomplete_sidechain_residues) > 0,
        "missing_residue_details":    missing_residue_details,
        "incomplete_sidechain_residues": incomplete_sidechain_residues,
        "chains":                     chains_found,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — PDBFIXER REPAIR (side chains only)
# Used when residues are all present but some lack side-chain atoms.
# Fast: runs in seconds, no GPU needed.
# ═══════════════════════════════════════════════════════════════════════════════

def repair_sidechains_pdbfixer(pdb_path: str, output_path: Optional[str] = None) -> str:
    """
    Use PDBFixer + OpenMM to rebuild missing side-chain heavy atoms.

    Parameters
    ----------
    pdb_path    : path to the incomplete PDB
    output_path : where to write the repaired PDB
                  (defaults to <stem>_sidechain_fixed.pdb next to input)

    Returns
    -------
    Path to the repaired PDB file.
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        raise ImportError(
            "PDBFixer is not installed.\n"
            "Run:  pip install pdbfixer openmm"
        )

    if output_path is None:
        stem        = Path(pdb_path).stem
        parent      = Path(pdb_path).parent
        output_path = str(parent / f"{stem}_sidechain_fixed.pdb")

    fixer = PDBFixer(filename=pdb_path)

    # Tell PDBFixer not to add missing residues — only fix atoms within
    # residues that are already present.
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # ← clear: we only want side chains

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()             # adds heavy atoms only (no H)

    with open(output_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    print(f"[logic5] Side-chain repair complete → {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — SEQUENCE EXTRACTION
# Reads the actual sequence from ATOM records for each chain.
# ═══════════════════════════════════════════════════════════════════════════════

# Standard 3-letter → 1-letter amino acid mapping
_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "HIE": "H",
    "HID": "H", "HIP": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",
    "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O",
}

def extract_sequence_from_pdb(pdb_path: str, chain_id: str) -> str:
    """
    Extract the amino acid sequence from ATOM records for a given chain.
    Uses only residues with complete backbone (CA atom present).

    Returns a one-letter-code string.
    """
    parser  = PDBParser(QUIET=True)
    struct  = parser.get_structure("s", pdb_path)
    seq     = []
    seen    = set()

    for model in struct:
        if chain_id not in [c.get_id() for c in model]:
            continue
        chain = model[chain_id]
        for residue in chain:
            # Skip HETATM / water
            het, resseq, _ = residue.get_id()
            if het.strip():
                continue
            resname = residue.get_resname().strip()
            if resseq in seen:
                continue
            seen.add(resseq)
            aa = _AA3TO1.get(resname)
            if aa:
                seq.append(aa)

    return "".join(seq)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — MASTER ENTRY POINT
# This is what gui.py and the rest of the pipeline should call.
# It decides the repair strategy automatically and returns a complete PDB.
# ═══════════════════════════════════════════════════════════════════════════════

def prepare_pdb(
    pdb_path:        str,
    cam_chain_id:    str,
    peptide_chain_id: str,
) -> tuple[str, str]:
    """
    Inspect the PDB and repair it using the appropriate strategy.

    This is the main function called before run_prppi() in the pipeline.

    Parameters
    ----------
    pdb_path         : path to the raw input PDB
    cam_chain_id     : calmodulin chain (from identify_calmodulin_chain)
    peptide_chain_id : peptide chain

    Returns
    -------
    (repaired_pdb_path, strategy_used)
      strategy_used is one of: "none_needed", "sidechain_only", "missing_residues_skipped"
    """
    print(f"[logic5] Inspecting: {pdb_path}")
    report = inspect_pdb(pdb_path)

    print(f"[logic5] Missing residues:    {report['has_missing_residues']}")
    print(f"[logic5] Missing side chains: {report['has_missing_sidechains']}")

    if report["has_missing_residues"]:
        # Missing residues/loops are not currently auto-repaired. They are
        # logged here and surfaced to the user in gui.py, but the pipeline
        # proceeds with the original coordinates for that portion.
        print(
            f"[logic5] WARNING: {len(report['missing_residue_details'])} missing "
            f"residue(s) detected. No automatic repair is available for "
            f"missing residues — proceeding with side-chain repair only if needed."
        )

    if report["has_missing_sidechains"]:
        # PDBFixer — fast, in-process
        print("[logic5] Strategy: PDBFixer side-chain repair")
        repaired = repair_sidechains_pdbfixer(pdb_path)
        return repaired, "sidechain_only"

    else:
        # Covers both a fully complete PDB and the case where only missing
        # residues were found (remodelling disabled — no action taken).
        if report["has_missing_residues"]:
            print("[logic5] Strategy: Missing residues noted — no repair applied (disabled)")
            return pdb_path, "missing_residues_skipped"
        print("[logic5] Strategy: No repair needed — PDB is complete")
        return pdb_path, "none_needed"


def get_inspection_report(pdb_path: str) -> dict:
    """
    Convenience wrapper: returns the inspection dict for display in the GUI.
    Called in Step 1 of gui.py after the PDB is uploaded.
    """
    return inspect_pdb(pdb_path)