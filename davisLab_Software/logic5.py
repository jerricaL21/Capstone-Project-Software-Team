"""
logic5.py — PDB Structural Repair
==================================
Sits between Step 1 (Data Input) and the prppi call in logic.py.

Responsibilities:
  1. Inspect an incoming PDB for missing side chains and/or missing residues.
  2. Rebuild missing internal residues and missing side-chain atoms via
     PDBFixer (fast, in-process).
  3. Pick the least-clashing rotamer for every rebuilt / completed residue
     using PyMOL's mutagenesis wizard.
  4. Return a clean, complete PDB path ready for fix_pdb_protonation()
     and run_prppi() in logic.py.

Notes
-----
* Chain-end (N/C-terminal) missing residues are  rebuilt by 
  (skip_termini=False) 
* Rebuilt loop coordinates are approximate (PDBFixer places them roughly
  along a straight line between the gap ends). Rotamer selection only
  resolves side-chain clashes, not backbone placement.

Dependencies
------------
  pip install pdbfixer openmm
  PyMOL (pymol2 or pymol) — already required by logic6.py / Step 5
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
      has_missing_sidechains (bool) — residues present but lacking side-chain atoms
      missing_residue_details (list of str) — human-readable descriptions
      incomplete_sidechain_residues (list of str)
      incomplete_sidechain_keys (list of (chain_id, resnum)) — structured form
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
        # GLY has no side chain — skipped in the loop below
        # SEC (selenocysteine) and PYL (pyrrolysine) are rare; omitted
    }

    incomplete_sidechain_residues = []
    incomplete_sidechain_keys = []
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
                    incomplete_sidechain_keys.append((cid, resid))

    return {
        "has_missing_residues":          len(missing_residue_details) > 0,
        "has_missing_sidechains":        len(incomplete_sidechain_residues) > 0,
        "missing_residue_details":       missing_residue_details,
        "incomplete_sidechain_residues": incomplete_sidechain_residues,
        "incomplete_sidechain_keys":     incomplete_sidechain_keys,
        "chains":                        chains_found,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — PDBFIXER REPAIR (missing residues + side chains)
# Fast: runs in seconds, no GPU needed.
# ═══════════════════════════════════════════════════════════════════════════════

def _residue_keys(pdb_path: str) -> set:
    """Return {(chain_id, resnum), ...} for all standard (non-hetero) residues."""
    s = PDBParser(QUIET=True).get_structure("s", pdb_path)
    return {(c.id, r.id[1]) for c in s[0] for r in c if r.id[0].strip() == ""}


def repair_structure_pdbfixer(
    pdb_path: str,
    output_path: Optional[str] = None,
    rebuild_missing: bool = True,
    skip_termini: bool = True,
):
    """
    Use PDBFixer + OpenMM to rebuild missing residues (internal gaps) and
    missing side-chain heavy atoms, then add hydrogens at pH 7.0.

    Parameters
    ----------
    pdb_path        : path to the incomplete PDB
    output_path     : where to write the repaired PDB
                      (defaults to <stem>_sidechain_fixed.pdb next to input)
    rebuild_missing : if False, only side-chain atoms are rebuilt
    skip_termini    : if True, missing residues at chain ends are not rebuilt

    Returns
    -------
    (output_path, rebuilt_keys)
      rebuilt_keys = sorted [(chain_id, resnum), ...] of residues that exist
      in the output but not in the input (i.e. newly built residues).
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
        p = Path(pdb_path)
        output_path = str(p.parent / f"{p.stem}_sidechain_fixed.pdb")

    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()

    if not rebuild_missing:
        fixer.missingResidues = {}
    elif skip_termini:
        chains = list(fixer.topology.chains())
        fixer.missingResidues = {
            (ci, ri): names
            for (ci, ri), names in fixer.missingResidues.items()
            if 0 < ri < len(list(chains[ci].residues()))   # internal gaps only
        }

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()           # builds missing residues + missing side-chain atoms
    fixer.addMissingHydrogens(7.0)

    with open(output_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    # Diffing residue IDs before/after is more robust than predicting
    # PDBFixer's numbering for the new residues.
    new_keys = _residue_keys(output_path) - _residue_keys(pdb_path)

    print(f"[logic5] Structure repair complete → {output_path} "
          f"({len(new_keys)} residue(s) newly built)")
    return output_path, sorted(new_keys)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — PYMOL ROTAMER SELECTION
# For each rebuilt / completed residue, let PyMOL's mutagenesis wizard
# generate its rotamer library ("mutating" the residue to itself) and keep the
# state with the lowest bump score (fewest steric clashes).
# ═══════════════════════════════════════════════════════════════════════════════

_HIS_LIKE = {"HIE": "HIS", "HID": "HIS", "HIP": "HIS"}


def optimize_rotamers_pymol(pdb_path: str, residues: list, output_path: Optional[str] = None):
    """
    Parameters
    ----------
    pdb_path    : PDB to optimize
    residues    : [(chain_id, resnum), ...]
    output_path : defaults to <stem>_rotamers.pdb next to input

    Returns
    -------
    (output_path, log_lines)
    """
    from logic6 import _start_pymol_session   # reuse the PyMOL session helper

    if output_path is None:
        p = Path(pdb_path)
        output_path = str(p.parent / f"{p.stem}_rotamers.pdb")

    log = []
    cmd, stop = _start_pymol_session()
    try:
        cmd.load(pdb_path, "target")
        for chain, resnum in residues:
            sel = f"target and chain {chain} and resi {resnum}"
            names = []
            cmd.iterate(sel + " and name CA", "names.append(resn)", space={"names": names})
            if not names:
                continue
            resn = _HIS_LIKE.get(names[0], names[0])
            if resn in ("GLY", "ALA"):          # no rotamers to choose between
                continue

            tmp = f"rot_{chain}{resnum}"
            try:
                cmd.select(tmp, sel)
                cmd.wizard("mutagenesis")
                cmd.refresh_wizard()
                wiz = cmd.get_wizard()
                wiz.set_mode(resn)
                wiz.do_select(tmp)

                # one score per rotamer state; lower = fewer clashes.
                # NOTE: bump_scores is an internal wizard attribute — verify
                # with print(dir(wiz)) on your PyMOL build.
                scores = list(getattr(wiz, "bump_scores", []))
                if scores:
                    best = min(range(len(scores)), key=scores.__getitem__) + 1
                    cmd.frame(best)
                    log.append(f"{chain}{resnum} {resn}: rotamer {best}/{len(scores)} "
                               f"(bump {scores[best - 1]:.2f})")
                else:
                    cmd.frame(1)
                    log.append(f"{chain}{resnum} {resn}: no scores, kept rotamer 1")
                wiz.apply()
                cmd.set_wizard()
            except Exception as e:
                log.append(f"{chain}{resnum}: failed ({e})")
                try:
                    cmd.set_wizard()
                except Exception:
                    pass
            finally:
                try:
                    cmd.delete(tmp)
                except Exception:
                    pass

        cmd.save(output_path, "target")
    finally:
        stop()

    for line in log:
        print(f"[logic5] {line}")
    return output_path, log


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — SEQUENCE EXTRACTION
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
# SECTION 5 — MASTER ENTRY POINT
# This is what gui.py and the rest of the pipeline should call.
# ═══════════════════════════════════════════════════════════════════════════════

def prepare_pdb(
    pdb_path:         str,
    cam_chain_id:     str,
    peptide_chain_id: str,
    return_details:   bool = False,
):
    """
    Inspect the PDB and repair it: PDBFixer rebuilds missing internal residues
    and side chains, then PyMOL picks the least-clashing rotamer for each
    rebuilt / completed residue.

    Returns
    -------
    return_details=False (default): (repaired_pdb_path, strategy_used)
    return_details=True           : (repaired_pdb_path, strategy_used, target_residues)

      strategy_used is one of:
        "none_needed"
        "sidechain_rotamer_optimized"
        "rebuilt_and_rotamer_optimized"
      target_residues = [(chain_id, resnum), ...] that were rebuilt/optimized
    """
    print(f"[logic5] Inspecting: {pdb_path}")
    report = inspect_pdb(pdb_path)

    print(f"[logic5] Missing residues:    {report['has_missing_residues']}")
    print(f"[logic5] Missing side chains: {report['has_missing_sidechains']}")

    needs_fix = report["has_missing_sidechains"] or report["has_missing_residues"]

    if not needs_fix:
        print("[logic5] Strategy: No repair needed — PDB is complete")
        result = (pdb_path, "none_needed", [])
    else:
        print("[logic5] Strategy: PDBFixer rebuild + PyMOL rotamer selection")
        repaired, rebuilt = repair_structure_pdbfixer(pdb_path, skip_termini=False)
        targets = sorted(set(rebuilt) | set(report["incomplete_sidechain_keys"]))
        final, _log = optimize_rotamers_pymol(repaired, targets)
        strategy = "rebuilt_and_rotamer_optimized" if rebuilt else "sidechain_rotamer_optimized"
        result = (final, strategy, targets)

    return result if return_details else result[:2]


def get_inspection_report(pdb_path: str) -> dict:
    """
    Convenience wrapper: returns the inspection dict for display in the GUI.
    Called in Step 1 of gui.py after the PDB is uploaded.
    """
    return inspect_pdb(pdb_path)