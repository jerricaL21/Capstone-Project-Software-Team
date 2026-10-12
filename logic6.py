"""
logic6.py — PyMOL-Based Structural Redesign
=============================================
Sits after Step 4 (Results & Analytics) in gui.py.
"""

import os
import re
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
import pandas as pd

MIN_SITES = 1

_SITE_RE = re.compile(r"^([A-Za-z]+)(\d+)=(\w+)$")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — PICK THE BEST MUTATION PER SCANNED SITE
# ═══════════════════════════════════════════════════════════════════════════════

def parse_residue_site(wt_full: str) -> Tuple[str, int, str]:
    """Parse a WT_Full residue-site identifier (e.g. "A11=GLU") into (chain, resnum, wt_aa3)."""
    m = _SITE_RE.match(str(wt_full).strip())
    if not m:
        raise ValueError(
            f"Could not parse residue site from '{wt_full}'. "
            f"Expected a format like 'A11=GLU'."
        )
    chain, resnum, wt_aa3 = m.groups()
    return chain, int(resnum), wt_aa3.upper()


def select_best_mutations(
    pivot_delta: pd.DataFrame, 
    selected_sites: Optional[List[str]] = None
) -> List[Dict[str, Any]]:
    """
    Select the best mutation per site from pivot_delta.
    If selected_sites is provided, filters for only those specific site identifiers.
    """
    if pivot_delta is None or pivot_delta.empty:
        raise ValueError("No mutation data available to select best mutations.")

    best_mutations = []
    for wt_full, row in pivot_delta.iterrows():
        site_str = str(wt_full).strip()
        
        # Filter by selected sites if provided
        if selected_sites is not None and site_str not in selected_sites:
            continue

        valid_scores = row.dropna()
        if valid_scores.empty:
            continue

        # Best mutation is the one with highest Delta K* score
        best_mut_aa = valid_scores.idxmax()
        best_delta = float(valid_scores[best_mut_aa])

        chain, resnum, wt_aa3 = parse_residue_site(site_str)

        best_mutations.append({
            "site": site_str,
            "chain": chain,
            "resnum": resnum,
            "wt_aa3": wt_aa3,
            "mutant_aa3": best_mut_aa,
            "delta": best_delta,
        })

    # Sort selected sites by highest Delta K* score descending
    best_mutations.sort(key=lambda x: x["delta"], reverse=True)

    return best_mutations


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — APPLY MUTATIONS IN PYMOL
# ═══════════════════════════════════════════════════════════════════════════════

def _start_pymol_session():
    """Return (cmd, stop_fn) using whichever PyMOL Python API is importable."""
    try:
        import pymol2
        session = pymol2.PyMOL()
        session.start()
        return session.cmd, session.stop
    except ImportError:
        pass

    try:
        import pymol
        from pymol import cmd as pymol_cmd
        if not getattr(pymol, "_logic6_launched", False):
            pymol.finish_launching(["pymol", "-qc"])  # -q quiet, -c no GUI
            pymol._logic6_launched = True
        return pymol_cmd, (lambda: None)
    except ImportError as e:
        import sys
        raise ImportError(
            "Could not import a PyMOL Python API ('pymol2' or 'pymol') in "
            f"this interpreter:\n  {sys.executable}\n\n"
            "This usually means PyMOL is installed in a DIFFERENT Python "
            "environment than the one running this app.\n"
            f"  • install PyMOL into this interpreter: {sys.executable} -m pip install pymol-open-source\n"
            "  • or launch this app from the environment where PyMOL is installed."
        ) from e


def apply_mutations_pymol(
    pdb_path: str,
    mutations: list[dict],
    output_path: Optional[str] = None,
):
    """
    Use PyMOL's mutagenesis wizard to introduce each mutation in `mutations`
    into `pdb_path`, then save the combined multi-site mutant structure.
    """
    if not mutations:
        raise ValueError("No mutations to apply.")

    if not pdb_path or not os.path.isfile(pdb_path):
        raise FileNotFoundError(
            f"Input PDB file not found at:\n  {pdb_path}\n"
        )

    if output_path is None:
        stem = Path(pdb_path).stem
        parent = Path(pdb_path).parent
        output_path = str(parent / f"{stem}_pymol_redesigned.pdb")

    obj_name = "target"
    applied, errors = [], []

    cmd, stop_session = _start_pymol_session()
    try:
        cmd.load(pdb_path, obj_name)

        loaded_objects = cmd.get_names("objects")
        if obj_name not in loaded_objects:
            raise RuntimeError(
                f"PyMOL did not create an object named '{obj_name}' after loading:\n  {pdb_path}"
            )

        atom_count = cmd.count_atoms(obj_name)
        if atom_count == 0:
            raise RuntimeError(f"PyMOL loaded '{obj_name}' from {pdb_path}, but it contains 0 atoms.")

        for i, m in enumerate(mutations):
            chain, resnum, mutant_aa3 = m["chain"], m["resnum"], m["mutant_aa3"]
            sel = f"{obj_name} and chain {chain} and resi {resnum}"
            tmp_sel_name = f"mut_sel_{i}"

            try:
                if obj_name not in cmd.get_names("objects"):
                    raise RuntimeError(f"Object '{obj_name}' no longer exists in PyMOL session.")

                if cmd.count_atoms(sel) == 0:
                    errors.append(f"Chain {chain} residue {resnum}: not found in structure — skipped.")
                    continue

                cmd.select(tmp_sel_name, sel)
                cmd.wizard("mutagenesis")
                cmd.refresh_wizard()
                cmd.get_wizard().set_mode(mutant_aa3)
                cmd.get_wizard().do_select(tmp_sel_name)
                cmd.frame(1)
                cmd.get_wizard().apply()
                cmd.set_wizard()
                applied.append(
                    f"{chain}{resnum} {m.get('wt_aa3', '?')} \u2192 {mutant_aa3} "
                    f"(\u0394K* = {m.get('delta', 0):+.3f})"
                )
            except Exception as e:
                errors.append(f"Chain {chain} residue {resnum} \u2192 {mutant_aa3}: {e}")
                try:
                    cmd.set_wizard()
                except Exception:
                    pass
            finally:
                try:
                    cmd.delete(tmp_sel_name)
                except Exception:
                    pass

        cmd.save(output_path, obj_name)

    finally:
        stop_session()

    return output_path, applied, errors


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — ONE-CALL ENTRY POINT FOR gui.py's STEP 5
# ═══════════════════════════════════════════════════════════════════════════════

def generate_redesigned_structure(
    pdb_path: str,
    pivot_delta: pd.DataFrame,
    selected_sites: Optional[List[str]] = None,
    output_path: Optional[str] = None,
) -> Tuple[str, List[Dict[str, Any]], List[str], List[str]]:
    """
    Entry point called by gui.py (Step 5).

    1. Selects best mutations for the user-selected sites from pivot_delta.
    2. Applies all of them to pdb_path in PyMOL via apply_mutations_pymol().
    3. Returns (output_path, mutations, applied, errors).
    """
    mutations = select_best_mutations(pivot_delta, selected_sites=selected_sites)
    output_path, applied, errors = apply_mutations_pymol(pdb_path, mutations, output_path)
    return output_path, mutations, applied, errors