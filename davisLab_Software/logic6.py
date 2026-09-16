"""
logic6.py — PyMOL-Based Structural Redesign
=============================================
Sits after Step 4 (Results & Analytics) in gui.py.

Workflow position
------------------
By the time this runs, the user has completed at least three single-site
mutation scans in Step 3 (the original scan plus 2 additional ones, as
requested each producing its own bbkstar_results_*.tsv under the OSPREY
output tree that logic4.analyze_results() already walks and aggregates).

Responsibilities
-----------------
  1. From the aggregated ΔK* results (logic4's pivot_delta table), pick
     the single best-scoring mutation for every residue site that was
     scanned.
  2. Drive PyMOL's mutagenesis wizard to introduce all of those
     substitutions into the working PDB structure, one residue at a time.
  3. Save and return the resulting multi-site "best mutant" PDB.

Dependencies
------------
  pip install pymol-open-source     (provides the `pymol2` module used here)
"""

import os
import re
from pathlib import Path
from typing import Optional

# Original single-site scan (Step 3, first pass) + 2 more requested = 3.
MIN_SITES = 3

# logic4.read_file_data() builds WT_Full / Mutation strings like "A11=GLU"
# (chain letters, residue number, "=", 3-letter amino acid code).
_SITE_RE = re.compile(r"^([A-Za-z]+)(\d+)=(\w+)$")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — PICK THE BEST MUTATION PER SCANNED SITE
# ═══════════════════════════════════════════════════════════════════════════════

def parse_residue_site(wt_full: str) -> tuple[str, int, str]:
    """
    Parse a WT_Full residue-site identifier (e.g. "A11=GLU") produced by
    logic4.read_file_data() into (chain, resnum, wt_aa3).
    """
    m = _SITE_RE.match(str(wt_full).strip())
    if not m:
        raise ValueError(
            f"Could not parse residue site from '{wt_full}'. "
            f"Expected a format like 'A11=GLU'."
        )
    chain, resnum, wt_aa3 = m.groups()
    return chain, int(resnum), wt_aa3.upper()


def select_best_mutations(pivot_delta, min_sites: int = MIN_SITES) -> list[dict]:
    """
    Given the ΔK* pivot table from logic4.analyze_results()
    (index=WT_Full residue site, columns=Mutation, values=Delta, higher
    delta = greater predicted binding-affinity improvement — same
    convention gui.py Step 4 uses for "Top ΔK* Improvement"), pick the
    single best-scoring mutation for every scanned site.

    Raises ValueError if fewer than `min_sites` distinct residue sites
    have data yet — the caller should send the user back to Step 3 to
    scan more residues before offering the PyMOL redesign.

    Returns a list of dicts, sorted best-first:
        {"site": "A11=GLU", "chain": "A", "resnum": 11,
         "wt_aa3": "GLU", "mutant_aa3": "ALA", "delta": 1.234}
    """
    if pivot_delta is None or pivot_delta.empty:
        raise ValueError("No mutation-scan data available yet.")

    n_sites = pivot_delta.shape[0]
    if n_sites < min_sites:
        remaining = min_sites - n_sites
        raise ValueError(
            f"Only {n_sites} residue site(s) have been scanned so far. "
            f"Please return to Step 3 and run OSPREY on {remaining} more "
            f"residue site(s) before generating the PyMOL redesign."
        )

    best_mutations = []
    seen_sites = {}  # (chain, resnum) -> index into best_mutations, for dedup
    for site in pivot_delta.index:
        row = pivot_delta.loc[site].dropna()
        if row.empty:
            continue

        best_mutant_aa3 = row.idxmax()
        best_delta = float(row.loc[best_mutant_aa3])
        chain, resnum, wt_aa3 = parse_residue_site(site)

        candidate = {
            "site":       site,
            "chain":      chain,
            "resnum":     resnum,
            "wt_aa3":     wt_aa3,
            "mutant_aa3": str(best_mutant_aa3).upper(),
            "delta":      best_delta,
        }

        # A residue site can appear more than once if it was scanned in more
        # than one Step 3 run. PyMOL's mutagenesis wizard doesn't reliably
        # support mutating the same residue twice in one session, so keep
        # only the best-scoring entry per (chain, resnum) instead of
        # attempting the mutation multiple times.
        key = (chain, resnum)
        if key in seen_sites:
            existing_idx = seen_sites[key]
            if candidate["delta"] > best_mutations[existing_idx]["delta"]:
                best_mutations[existing_idx] = candidate
            continue

        seen_sites[key] = len(best_mutations)
        best_mutations.append(candidate)

    best_mutations.sort(key=lambda d: d["delta"], reverse=True)
    return best_mutations


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — APPLY MUTATIONS IN PYMOL
# ═══════════════════════════════════════════════════════════════════════════════

def _start_pymol_session():
    """
    Return (cmd, stop_fn) using whichever PyMOL Python API is importable
    in THIS interpreter.

    Tries, in order:
      1. `pymol2` — the modern, sandboxed multi-instance API
         (`pip install pymol-open-source` and most conda-forge builds
         provide this).
      2. Classic `pymol` module with `pymol.finish_launching(["-qc"])` —
         older builds, or some packaging of `pymol-open-source`, only
         ship this API.

    If neither import succeeds, PyMOL is installed somewhere this
    interpreter can't see (a different conda env / venv, or a standalone
    PyMOL app not on this Python's sys.path) rather than truly absent,
    so the error message points at the mismatch instead of assuming
    PyMOL needs to be installed from scratch.
    """
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
            "environment than the one running this app (e.g. a conda env, "
            "or a standalone PyMOL app not on this interpreter's path). "
            "Either:\n"
            f"  • install PyMOL into this interpreter: {sys.executable} -m pip install pymol-open-source\n"
            "  • or launch this app from the environment where PyMOL is "
            "already installed (e.g. `conda activate <env>` first)."
        ) from e


def apply_mutations_pymol(
    pdb_path: str,
    mutations: list[dict],
    output_path: Optional[str] = None,
):
    """
    Use PyMOL's mutagenesis wizard to introduce each mutation in
    `mutations` (as produced by select_best_mutations()) into `pdb_path`,
    then save the combined multi-site mutant structure.

    Each site is mutated to its top-ranked rotamer (wizard frame 1) —
    the same default PyMOL applies interactively if you don't hand-pick
    an alternate rotamer.

    Returns (output_path, applied, errors):
        applied — list of human-readable strings, one per mutation made
        errors  — list of human-readable strings for any site that could
                  not be mutated (e.g. chain/residue not found)
    """
    if not mutations:
        raise ValueError("No mutations to apply.")

    # -----------------------------------------------------------------
    # Fail loudly and early if the input PDB doesn't actually exist.
    # Without this check, cmd.load() can fail silently in PyMOL and the
    # first visible symptom is a confusing "Invalid selection name"
    # error several lines later, with no indication of the real cause.
    # -----------------------------------------------------------------
    if not pdb_path or not os.path.isfile(pdb_path):
        raise FileNotFoundError(
            f"Input PDB file not found at:\n  {pdb_path}\n\n"
            "Check that this path is correct and the file still exists "
            "(it may have been generated in a temp/output folder that "
            "was moved, cleared, or never created)."
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

        # -----------------------------------------------------------------
        # Verify the load actually worked. PyMOL's cmd.load() does not
        # always raise an exception on failure (e.g. malformed PDB,
        # permissions issue, or unexpected file content) - it can just
        # print a warning and skip creating the object. If we don't check
        # here, every later command referencing obj_name fails instead
        # with a confusing "Invalid selection name" error that hides the
        # real cause.
        # -----------------------------------------------------------------
        loaded_objects = cmd.get_names("objects")
        if obj_name not in loaded_objects:
            raise RuntimeError(
                f"PyMOL did not create an object named '{obj_name}' after "
                f"loading:\n  {pdb_path}\n\n"
                f"Objects currently loaded in this PyMOL session: {loaded_objects}\n"
                "This usually means the file failed to load (e.g. it's "
                "empty, malformed, or not actually a valid PDB file)."
            )

        atom_count = cmd.count_atoms(obj_name)
        if atom_count == 0:
            raise RuntimeError(
                f"PyMOL loaded '{obj_name}' from:\n  {pdb_path}\n"
                "but it contains 0 atoms. The file may be empty or corrupted."
            )

        for i, m in enumerate(mutations):
            chain, resnum, mutant_aa3 = m["chain"], m["resnum"], m["mutant_aa3"]
            sel = f"{obj_name} and chain {chain} and resi {resnum}"

            # -----------------------------------------------------------------
            # IMPORTANT: never hand a compound "and"-based expression like
            # `sel` directly to the mutagenesis wizard's do_select(). Internally,
            # do_select() ends by calling cmd.delete() on whatever string it was
            # given, and cmd.delete() treats its argument as a whitespace-
            # separated list of NAMES to delete — not a selection expression.
            # Since `sel` contains the literal token "target" (our object name),
            # that call silently deletes the whole working structure as a side
            # effect of merely selecting a residue.
            #
            # Fix: pre-build a small, uniquely-named, single-token selection
            # with cmd.select(), and pass ONLY that name to do_select(). Then
            # whatever internal cleanup happens only removes the temporary
            # selection, never the object itself.
            # -----------------------------------------------------------------
            tmp_sel_name = f"mut_sel_{i}"

            try:
                # Confirm the base object is still intact before each
                # mutation attempt. In rare cases PyMOL's mutagenesis
                # wizard can leave the object in an unexpected state after
                # a prior mutation; catching that here turns it into a
                # per-mutation error instead of crashing the whole run.
                if obj_name not in cmd.get_names("objects"):
                    raise RuntimeError(
                        f"Object '{obj_name}' no longer exists in this "
                        f"PyMOL session (it may have been altered by a "
                        f"previous mutation step)."
                    )

                if cmd.count_atoms(sel) == 0:
                    errors.append(
                        f"Chain {chain} residue {resnum}: not found in the "
                        f"structure — mutation skipped."
                    )
                    continue

                # Build the temp named selection from the compound expression
                # ONCE here, on the real object. This is safe: cmd.select()
                # takes (name, selection-expression) as separate arguments,
                # so `sel`'s "and"-tokens are parsed correctly as a selection,
                # not as a list of names to delete.
                cmd.select(tmp_sel_name, sel)

                cmd.wizard("mutagenesis")
                cmd.refresh_wizard()
                cmd.get_wizard().set_mode(mutant_aa3)
                # Pass only the single-token temp selection name — never the
                # raw compound `sel` string — so that do_select()'s internal
                # cmd.delete() cleanup can only remove tmp_sel_name, not the
                # "target" object.
                cmd.get_wizard().do_select(tmp_sel_name)
                cmd.frame(1)               # top-ranked rotamer
                cmd.get_wizard().apply()
                cmd.set_wizard()           # clear the wizard before the next site
                applied.append(
                    f"{chain}{resnum} {m.get('wt_aa3', '?')} \u2192 {mutant_aa3} "
                    f"(\u0394K* = {m.get('delta', 0):+.3f})"
                )
            except Exception as e:
                errors.append(
                    f"Chain {chain} residue {resnum} \u2192 {mutant_aa3}: {e}"
                )
                try:
                    cmd.set_wizard()
                except Exception:
                    pass
            finally:
                # do_select()'s internal cleanup normally deletes
                # tmp_sel_name already, but that cleanup is exactly the
                # behavior we can't fully trust here — so also try to
                # remove it ourselves. It's just a named selection (never
                # the "target" object), so deleting it again if it's
                # already gone is harmless; ignore errors either way.
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
    pivot_delta,
    output_path: Optional[str] = None,
    min_sites: int = MIN_SITES,
):
    """
    THIS IS THE FUNCTION CALLED BY gui.py (Step 5).

    1. Picks the best mutation per scanned residue site from pivot_delta.
    2. Applies all of them to pdb_path in PyMOL.
    3. Returns (output_path, mutations, applied, errors).

    Raises ValueError if fewer than `min_sites` residue sites have been
    scanned yet (see select_best_mutations for the exact message shown
    to the user).
    """
    mutations = select_best_mutations(pivot_delta, min_sites=min_sites)
    output_path, applied, errors = apply_mutations_pymol(pdb_path, mutations, output_path)
    return output_path, mutations, applied, errors