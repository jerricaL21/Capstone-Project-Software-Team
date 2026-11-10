import streamlit as st
import pandas as pd
import time
import os
import importlib
import shutil
from install_pymol import ensure_pymol

# ── Ensure PyMOL is available in this interpreter (cached so it only
#    checks/installs once per session, not on every Streamlit rerun) ──
@st.cache_resource
def _pymol_ready():
    return ensure_pymol()

pymol_ok = _pymol_ready()
if not pymol_ok:
    st.sidebar.error(
        "⚠️ PyMOL not available in this Python — PDB repair (Step 1) will fail for "
        "structures with missing atoms/residues. Launch the app with run_app.bat "
        "(uses the pymol-env conda environment)."
    )

# Importing the script files
import logic
import logic2
import logic3
import logic4
import logic5

importlib.reload(logic)
importlib.reload(logic2)
importlib.reload(logic3)
importlib.reload(logic4)
importlib.reload(logic5)

# Used in Step 3 to figure out which three-letter amino acid code corresponds to wild-type residue code
AA_1TO3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}

# Page Configuration
st.set_page_config(page_title="Davis Lab Redesign", layout="wide")

# UI Header
st.title("CaM-RyR2 Protein Redesign Portal")
st.markdown("""
**Status:** Unified Workflow v1.0 (Draft)  
*Meeting ISO/IEC 25010:2023 Efficiency Standards*
""")

st.sidebar.header("Pipeline Controls")
ABOUT = "About"
workflow_step = st.sidebar.radio(
    "Navigate Workflow", 
    [ABOUT, "1. Data Input", "2. JSON Processing", "3. OSPREY Execution", "4. Results & Analytics"],
    key="workflow_step"
)

# ── "Next step" button shown at the bottom of each page ──────────────────────
def _go_to_step(target):
    st.session_state["workflow_step"] = target

def next_step_button(target, label, ready=True):
    if not ready:
        return
    st.markdown("---")
    st.button(f"{label}  ➡️", key=f"next_to_{target}", on_click=_go_to_step,
              args=(target,), type="primary", use_container_width=True)

# ═══════════════════════════════════════════════════════════════════════════════
# Step 0: About Section
# ═══════════════════════════════════════════════════════════════════════════════
if workflow_step == ABOUT:
    st.header("About the Protein Redesign Portal")

    st.markdown("""
    ### Overview

    The **Protein Redesign Portal** is a Streamlit-based application developed to simplify the process of evaluating protein-protein interactions and screening mutations. 
    It brings together several computational tools into one interface, allowing users to upload a protein structure, identify interacting residues, prepare input files, 
    run OSPREY calculations, and analyze mutation scores.

    ---

    ### Workflow Overview

    The portal guides users through four main stages:

    **1. Data Input & Structure Preparation**
    - Upload a PDB file and identify the protein chains involved in the complex.
    - Manually select the chains when automatic identification is ambiguous or multiple chains are present.
    - Check and repair missing residues or incomplete side chains using PDBFixer and PyMOL.
    - Run PRPPI to identify interacting residues between protein chains and generate a JSON file.
        - **JSON:** A structured file format used to store interaction data for further processing.

    **2. JSON Processing**
    - Process the JSON file
    - Define the residue ranges used for protein design.

    **3. OSPREY Execution**
    - Select mutation sites and amino acid substitutions.
    - Configure and run OSPREY BBK* calculations.
        - **OSPREY:** A protein design software package that evaluates how mutations may affect protein-protein binding.    

    **4. Results & Analysis**
    - Analyze and compare mutant and wild-type K* scores using ΔK*.
    - View results through heatmaps and summary tables.
    - Export results as CSV files.

    ---

    ### Technology Stack

    | Component | Technology |
    |---|---|
    | User Interface | Streamlit |
    | Structural Analysis | Biopython |
    | Structure Preparation | PDBFixer, OpenMM |
    | Protein Design | OSPREY (BBK*) |
    | Interaction Analysis | PRPPI |
    | Data Processing | Pandas, NumPy |

    """)

    next_step_button("1. Data Input", "Begin with Data Input")

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Data Input
# ═══════════════════════════════════════════════════════════════════════════════
elif workflow_step == "1. Data Input":
    st.header("Step 1: Structural Data Input")

    uploaded_pdb = st.file_uploader("Upload Target PDB File", type=['pdb'], key="pdb_uploader")

    # ── A different PDB was uploaded: reset everything, then load the new file ──
    if uploaded_pdb is not None and uploaded_pdb.name != st.session_state.get("pdb_name"):
        _PROTECTED = {"workflow_step", "pdb_uploader"}
        for k in list(st.session_state.keys()):
            if k not in _PROTECTED:
                del st.session_state[k]

        os.makedirs("GUI_results", exist_ok=True)
        # Create a subfolder inside GUI_results named after the PDB file (without .pdb)
        pdb_stem = os.path.splitext(uploaded_pdb.name)[0]
        pdb_sub_dir = os.path.join("GUI_results", pdb_stem)
        
        # Creates the paths for the subfolders inside the PDB folder
        pdb_files_dir = os.path.join(pdb_sub_dir, "pdb_files")
        json_dir = os.path.join(pdb_sub_dir, "prppi_json")
        osprey_dir = os.path.join(pdb_sub_dir, "osprey_output")
        
        # Creates the actual subfolders inside the PDB folder
        os.makedirs(pdb_files_dir, exist_ok=True)
        os.makedirs(json_dir, exist_ok=True)
        os.makedirs(osprey_dir, exist_ok=True)
        
        # Save the original uploaded PDB inside the 'pdb_files' folder
        new_path = os.path.abspath(os.path.join(pdb_files_dir, uploaded_pdb.name))
        with open(new_path, "wb") as f:
            f.write(uploaded_pdb.getvalue())
        
        # Storing the location of each folder path
        
        # Save these paths so other parts of the app can find them
        st.session_state["pdb_files_dir"] = pdb_files_dir
        st.session_state["json_dir"]      = json_dir
        st.session_state["osprey_dir"]    = osprey_dir
        
        st.session_state["pdb_name"]  = uploaded_pdb.name
        st.session_state["orig_path"] = new_path   
        st.session_state["path"]      = new_path   
        st.session_state["pdb_text"]  = uploaded_pdb.getvalue().decode("utf-8", errors="ignore")

    pdb_name = st.session_state.get("pdb_name")

    if pdb_name is None:
        st.info("Upload a PDB file to begin.")
    else:
        path     = st.session_state["orig_path"]
        pdb_text = st.session_state["pdb_text"]

        if uploaded_pdb is None:
            st.caption(f"Currently loaded: **{pdb_name}**. Upload a different PDB to start over.")

        st.success(f"Successfully loaded: {pdb_name}")
        st.write(f"Saved file path: {path}")

        with st.expander("📄 View PDB File Header"):
            st.code(pdb_text[:500] + "...", language="text")

        # ── Chain identification: only runs once per uploaded file ────────
        if st.session_state.get("identified_for") != pdb_name:
            with st.spinner("Identifying chains..."):
                try:
                    chain_results = logic.identify_calmodulin_chain(path)
                    st.session_state["chain_results"]  = chain_results
                    st.session_state["identified_for"] = pdb_name
                except Exception as e:
                    st.error(f"Chain identification failed: {e}")
                    chain_results = {}
        else:
            chain_results = st.session_state["chain_results"]

        # ── Display chain identification results ──────────────────────────
        if chain_results:
            st.subheader("🔍 Chain Identification")

            with st.expander("View chain identification results", expanded=True):
                for chain_id, info in chain_results.items():
                    mol_name = info.get("molecule_name", f"Chain {chain_id}")
                    if info["is_calmodulin"]:
                        if info.get("identity") is not None:
                            identity_msg = f" ({info['identity']}% identity to human CaM, {info['length']} residues)"
                        else:
                            identity_msg = f" ({info['length']} residues)"
                        st.success(f"Chain **{chain_id}** ({mol_name}) → CALMODULIN{identity_msg}")
                    else:
                        st.info(f"Chain **{chain_id}** ({mol_name}) → Ligand / peptide ({info['length']} residues)")

            all_chains     = list(chain_results.keys())
            cam_candidates = [c for c in all_chains if chain_results[c]["is_calmodulin"]]
            other_chains   = [c for c in all_chains if not chain_results[c]["is_calmodulin"]]

            multi_chain  = len(all_chains) > 2
            needs_choice = multi_chain or len(cam_candidates) != 1

            # Only shown for PDB files with more than 2 chains
            if multi_chain:
                st.info("Please review your PDB file to understand its contents before selecting the chains to submit to PRPPI.")

            def _chain_label(c):
                i = chain_results[c]
                tag = "CaM" if i["is_calmodulin"] else "other"
                ident = f", {i['identity']}% identity" if i.get("identity") is not None else ""
                return f"{c} — {tag} ({i['length']} residues{ident})"

            if len(all_chains) < 2:
                st.error("At least 2 chains are required (calmodulin + partner).")
                cam_chain, ligand_chain = None, None
            elif needs_choice:
                st.warning(
                    "Multiple chains detected (or CaM could not be resolved automatically). "
                    "Please select which chain is the **main protein chain** and which is the "
                    "**ligand/partner chain**."
                )

                def _rank(c):
                    ident = chain_results[c].get("identity")
                    return (ident if ident is not None else -1, chain_results[c]["length"])

                default_cam = max(cam_candidates, key=_rank) if cam_candidates else all_chains[0]

                # Previously saved choices (these survive leaving and returning to this page)
                cam_prev = st.session_state.get("cam_chain")
                lig_prev = st.session_state.get("ligand_chain")

                if cam_prev in all_chains:
                    cam_index = all_chains.index(cam_prev)
                elif multi_chain:
                    cam_index = None          # more than 2 chains: selection required
                else:
                    cam_index = all_chains.index(default_cam)

                col_a, col_b = st.columns(2)
                with col_a:
                    cam_chain = st.selectbox(
                        "Main protein chain",
                        all_chains,
                        index=cam_index,
                        placeholder="Select calmodulin chain",
                        format_func=_chain_label,
                    )
                with col_b:
                    lig_options = [c for c in all_chains if c != cam_chain]
                    if lig_prev in lig_options:
                        lig_index = lig_options.index(lig_prev)
                    elif multi_chain:
                        lig_index = None      # more than 2 chains: selection required
                    else:
                        default_lig = next((c for c in lig_options if c in other_chains), lig_options[0])
                        lig_index = lig_options.index(default_lig)
                    ligand_chain = st.selectbox(
                        "Ligand / partner chain",
                        lig_options,
                        index=lig_index,
                        placeholder="Select peptide / ligand chain",
                        format_func=_chain_label,
                    )
            else:
                cam_chain, ligand_chain = cam_candidates[0], other_chains[0]

            # Save so the Run Analysis button (next rerun) and later steps can read them
            st.session_state["cam_chain"]    = cam_chain
            st.session_state["ligand_chain"] = ligand_chain

        # ── Gate: for >2 chains, wait for a selection before inspecting ───
        if chain_results and len(chain_results) > 2 and not (cam_chain and ligand_chain):
            st.info("Select the calmodulin and peptide/ligand chains above to continue.")
            st.stop()      

        # ── PDB Inspection ────────────────────────────────────────────────
        # Runs logic5.get_inspection_report() and shows the user what problems were found before the repair step.
        if st.session_state.get("inspected_for") != pdb_name:
            with st.spinner("Inspecting PDB for structural completeness..."):
                try:
                    report = logic5.get_inspection_report(path)
                    st.session_state["inspection_report"] = report
                    st.session_state["inspected_for"]     = pdb_name
                except Exception as e:
                    st.error(f"PDB inspection failed: {e}")
                    report = None
        else:
            report = st.session_state.get("inspection_report")

        if report:
            st.subheader("🔬 Structural Completeness Check")

            col1, col2 = st.columns(2)
            with col1:
                if report["has_missing_residues"]:
                    st.error(
                        f"❌ Missing residues detected "
                        f"({len(report['missing_residue_details'])} residues absent from coordinates)"
                    )
                    with st.expander("View missing residues"):
                        for r in report["missing_residue_details"]:
                            st.text(r)
                else:
                    st.success("✅ No missing residues")

            with col2:
                if report["has_missing_sidechains"]:
                    st.warning(
                        f"⚠️ Incomplete side chains detected "
                        f"({len(report['incomplete_sidechain_residues'])} residues affected)"
                    )
                    with st.expander("View incomplete side chains"):
                        for r in report["incomplete_sidechain_residues"]:
                            st.text(r)
                else:
                    st.success("✅ All side chains complete")

            # Describe what repair strategy will be used
            if report["has_missing_residues"] or report["has_missing_sidechains"]:
                msg = "🔧 **Repair strategy: PDBFixer + PyMOL rotamer selection.** "
                if report["has_missing_residues"]:
                    msg += ("Missing internal residues will be rebuilt. "
                            "Rebuilt loop coordinates are approximate, so treat contacts there with caution. ")
                if report["has_missing_sidechains"]:
                    msg += "Incomplete side chains will be completed. "
                msg += "PyMOL then picks the lowest-clash rotamer for each rebuilt residue."
                st.warning(msg)
            else:
                st.success("✅ PDB is complete — no repair needed before analysis.")

        # ── Run Analysis button ───────────────────────────────────────────
        if st.button("▶ Run Analysis"):

            cam_chain    = st.session_state.get("cam_chain")
            ligand_chain = st.session_state.get("ligand_chain")
            repair_info  = None

            # Step A: PDB repair (runs logic5 before prppi)
            if report and (report["has_missing_residues"] or report["has_missing_sidechains"]):
                if not cam_chain or not ligand_chain:
                    st.error(
                        "❌ Chain identification is required before repair. "
                        "Please ensure calmodulin and peptide chains were identified above."
                    )
                    st.stop()

                with st.spinner("Repairing PDB structure (PDBFixer + PyMOL)..."):
                    try:
                        repaired_path, strategy, rebuilt = logic5.prepare_pdb(
                            pdb_path         = path,
                            cam_chain_id     = cam_chain,
                            peptide_chain_id = ligand_chain,
                            return_details   = True,
                        )
                        st.session_state["path"] = repaired_path
                        repair_info = {"strategy": strategy, "repaired_path": repaired_path, "rebuilt": rebuilt}
                        path = repaired_path

                    except Exception as e:
                        st.error(f"❌ PDB repair failed: {e}")
                        st.stop()

            # Step B: Fix protonation (existing logic.py step)
            with st.spinner("Fixing PDB protonation..."):
                fixed_pdb = logic.fix_pdb_protonation(path)
                st.session_state["path"] = fixed_pdb

            # Step C: Run prppi (existing logic.py step)
            with st.spinner("Running PRPPI analysis..."):
                json_path = logic.run_prppi(
                    fixed_pdb, cutoff=5.0,
                    cam_info=st.session_state.get("chain_results"),
                    cam_chain=st.session_state.get("cam_chain"),
                    ligand_chain=st.session_state.get("ligand_chain"),
                )

            if json_path:
                # Copy JSON file into the prppi_json folder
                target_json_path = os.path.join(
                    st.session_state.get("json_dir", os.path.dirname(json_path)),
                    os.path.basename(json_path),
                )
                if os.path.abspath(json_path) != os.path.abspath(target_json_path):
                    shutil.copy(json_path, target_json_path)
                    json_path = target_json_path

                st.session_state["json_path"]       = json_path
                st.session_state["step1_done_for"]  = pdb_name
                st.session_state["analysis_result"] = {"repair": repair_info, "json_path": json_path}
                st.session_state.pop("step2_done", None)
                st.session_state.pop("step3_done", None)
                st.session_state.pop("step3_result", None)
                st.session_state.pop("osprey_scripts", None)
                st.session_state["scanned_residues"] = []

        # ── Show the latest analysis results (persists across steps) ──────
        res = st.session_state.get("analysis_result")
        if res:
            if res["repair"]:
                strategy_labels = {
                    "none_needed":                   "No repair was needed.",
                    "sidechain_rotamer_optimized":   "Side chains rebuilt with PDBFixer; rotamers chosen in PyMOL.",
                    "rebuilt_and_rotamer_optimized": "Missing residues and side chains rebuilt with PDBFixer; rotamers chosen in PyMOL.",
                }
                st.success(
                    f"✅ Repair complete: {strategy_labels.get(res['repair']['strategy'], res['repair']['strategy'])}\n\n"
                    f"Repaired file: `{res['repair']['repaired_path']}`"
                )
                if res["repair"]["rebuilt"]:
                    with st.expander(f"View {len(res['repair']['rebuilt'])} rebuilt / rotamer-optimized residues"):
                        st.text(", ".join(f"{c}{r}" for c, r in res["repair"]["rebuilt"]))
            else:
                st.info("PDB is complete — skipping repair step.")

            st.success(f"Analysis complete! Results saved to: {res['json_path']}")
            with st.expander("🔬 View JSON Contents"):
                import json as jsonlib
                with open(res["json_path"]) as f:
                    st.json(jsonlib.load(f))

    next_step_button("2. JSON Processing", "Click to head to JSON Processing",
                     ready=(pdb_name is not None and st.session_state.get("step1_done_for") == pdb_name))

# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: JSON Processing
# ═══════════════════════════════════════════════════════════════════════════════
elif workflow_step == "2. JSON Processing":
    st.header("Step 2: JSON Processing")

    json_path = st.session_state.get("json_path", None)

    if json_path:
        st.session_state["json_paths"] = [json_path]
        st.session_state["json_path"] = json_path
    else:
        uploaded_jsons = st.file_uploader(
            "Upload Interaction JSON(s)",
            type=['json'],
            accept_multiple_files=True
        )
        if uploaded_jsons:
            os.makedirs("uploaded_jsons", exist_ok=True)
            saved = []
            for uploaded_json in uploaded_jsons:
                save_path = os.path.abspath(os.path.join("uploaded_jsons", uploaded_json.name))
                with open(save_path, "wb") as f:
                    f.write(uploaded_json.getvalue())
                saved.append(save_path)
            st.session_state["json_paths"] = saved
            st.success(f"{len(saved)} JSON file(s) loaded:")
            for p in saved:
                st.write(f"• {p}")

    json_paths = st.session_state.get("json_paths", [])

    if json_paths:
        import logic2
        st.subheader("Set Residue Range Boundaries")

        cam_chain     = st.session_state.get("cam_chain", None)
        ligand_chain  = st.session_state.get("ligand_chain", None)
        chain_results = st.session_state.get("chain_results", {})

        pdb_path_for_residues = st.session_state.get("path", None)

        residue_map = {}
        if pdb_path_for_residues and os.path.exists(pdb_path_for_residues):
            try:
                residue_map = logic.get_residue_numbers(pdb_path_for_residues)
            except Exception:
                pass

        def get_residue_options(chain_id):
            if chain_id and residue_map:
                options = residue_map.get(chain_id, [])
                if options:
                    return options
            if chain_id and chain_id in chain_results:
                length = chain_results[chain_id]["length"]
                return [f"{chain_id}{i}" for i in range(1, length + 1)]
            return []

        cam_options    = get_residue_options(cam_chain)
        ligand_options = get_residue_options(ligand_chain)

        if cam_chain and ligand_chain:
            col1, col2 = st.columns(2)
            with col1:
                st.success(f"Calmodulin chain: **{cam_chain}** ({chain_results[cam_chain]['length']} residues)")
            with col2:
                st.info(f"Ligand chain: **{ligand_chain}** ({chain_results[ligand_chain]['length']} residues)")
        else:
            st.warning("⚠️ Chain identification not found. Please complete Step 1 first, or enter ranges manually.")

        def _idx(options, key, default):
            val = st.session_state.get(key)
            return options.index(val) if val in options else default

        with st.form("residue_form"):

            st.markdown("**Calmodulin (Mutant) Range**")
            if cam_options:
                col1, col2 = st.columns(2)
                with col1:
                    mutant_start = st.selectbox(
                        "Mutant Beginning",
                        options=cam_options,
                        index=_idx(cam_options, "saved_mutant_start", 0)
                    )
                with col2:
                    mutant_end = st.selectbox(
                        "Mutant Ending",
                        options=cam_options,
                        index=_idx(cam_options, "saved_mutant_end", len(cam_options) - 1)
                    )
            else:
                col1, col2 = st.columns(2)
                with col1:
                    mutant_start = st.text_input("Mutant Beginning", value=st.session_state.get("saved_mutant_start", ""))
                with col2:
                    mutant_end = st.text_input("Mutant Ending", value=st.session_state.get("saved_mutant_end", ""))

            st.markdown("**Ligand Range**")
            if ligand_options:
                col1, col2 = st.columns(2)
                with col1:
                    ligand_start = st.selectbox(
                        "Ligand Beginning",
                        options=ligand_options,
                        index=_idx(ligand_options, "saved_ligand_start", 0)
                    )
                with col2:
                    ligand_end = st.selectbox(
                        "Ligand Ending",
                        options=ligand_options,
                        index=_idx(ligand_options, "saved_ligand_end", len(ligand_options) - 1)
                    )
            else:
                col1, col2 = st.columns(2)
                with col1:
                    ligand_start = st.text_input("Ligand Beginning", value=st.session_state.get("saved_ligand_start", ""))
                with col2:
                    ligand_end = st.text_input("Ligand Ending", value=st.session_state.get("saved_ligand_end", ""))

            submitted = st.form_submit_button("▶ Process JSON(s)")

        if submitted:
            if not all([mutant_start, mutant_end, ligand_start, ligand_end]):
                st.error("❌ Please fill in all four residue range fields before proceeding.")
            else:
                for json_file in json_paths:
                    logic2.add_variables_to_json(json_file, mutant_start, mutant_end, ligand_start, ligand_end)
                st.session_state.update({
                    "saved_mutant_start": mutant_start, "saved_mutant_end": mutant_end,
                    "saved_ligand_start": ligand_start, "saved_ligand_end": ligand_end,
                    "step2_done": True,
                })
                st.success("All JSON files processed!")
        elif st.session_state.get("step2_done"):
            st.success("All JSON files processed!")

    next_step_button("3. OSPREY Execution", "Click to head to OSPREY Execution", ready=st.session_state.get("step2_done", False))

# ═══════════════════════════════════════════════════════════════════════════════
# Step 3: OSPREY Execution
# ═══════════════════════════════════════════════════════════════════════════════
elif workflow_step == "3. OSPREY Execution":
    st.header("Step 3: OSPREY Simulation Setup")
    
    json_path = st.session_state.get("json_path", None)
    pdb_path  = st.session_state.get("path", None)

    if not json_path or not pdb_path:
        st.warning("Please complete Steps 1 and 2 first.")
        st.stop()

    import json as jsonlib
    with open(json_path) as f:
        data = jsonlib.load(f)
    residues = [k for k in data.keys() if k != "Complex_Size"]
        
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Design Site")
        prev_res = st.session_state.get("step3_target_res")
        scanned = st.session_state.get("scanned_residues", [])
        target_res = st.selectbox(
            "Select Residue to Mutate", residues,
            index=residues.index(prev_res) if prev_res in residues else 0,
            format_func=lambda r: f"{r}  ✅ scanned" if r in scanned else r,
        )
        if scanned:
            st.caption(f"Already scanned this session: {', '.join(scanned)}. "
                       "Pick another residue to keep adding to the results.")

        # A different residue was chosen: clear the amino acid picks and the previous run's results
        if prev_res is not None and prev_res != target_res:
            st.session_state["saved_mutations"] = []
            st.session_state.pop("step3_result", None)
            st.session_state["step3_done"] = False
        st.session_state["step3_target_res"] = target_res

        all_aa_options = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'GLU', 'TYR', 'ASP', 'ARG', 'ASN', 'CYS', 'GLN', 'GLY', 'HIS', 'LYS', 'SER', 'THR']
        wt_letter = target_res[-1].upper() if target_res else None
        wt_three_letter = AA_1TO3.get(wt_letter)
        if wt_three_letter:
            mutation_options = [aa for aa in all_aa_options if aa != wt_three_letter]
        else:
            mutation_options = all_aa_options
        mutations = st.multiselect(
            "Select Amino Acids", mutation_options,
            default=[m for m in st.session_state.get("saved_mutations", []) if m in mutation_options],
        )

    with col2:
        st.subheader("Parameters")
        epsilon = st.slider("K* Precision (Epsilon)", 0.01, 1.0, st.session_state.get("saved_epsilon", 0.99))
        cores = st.number_input("CPU Cores", 1, 16, st.session_state.get("saved_cores", 4))
        gpus = st.number_input("GPUs", min_value=0, max_value=8, value=st.session_state.get("saved_gpus", 0), step=1)
        if gpus > 0:
            streams_per_gpu = st.number_input("Streams per GPU", min_value=0, max_value=256,
                                              value=st.session_state.get("saved_streams", 82), step=1)
            st.session_state["saved_streams"] = streams_per_gpu
        else:
            streams_per_gpu = 0
            st.info("Streams per GPU set to 0 (no GPU detected)")
        heap_size = st.number_input("Heap Size (MiB)", 1000, 500000, st.session_state.get("saved_heap", 100000))
        garbage_size = st.number_input("Garbage Size (MiB)", 512, 50000, st.session_state.get("saved_garbage", 8192))

    # Remember the current inputs so they are still here when you come back to this step
    st.session_state.update({
        "saved_mutations": mutations,
        "saved_epsilon":   epsilon,
        "saved_cores":     cores,
        "saved_gpus":      gpus,
        "saved_heap":      heap_size,
        "saved_garbage":   garbage_size,
    })

    if st.button("Generate & Run OSPREY Script"):
        st.session_state["step3_done"] = False 

        #-------At least one amino acid must be selected----#
        if not mutations:
            st.error("❌ Please select at least one amino acid before running OSPREY.")
            st.stop()
        try:
            with st.spinner("Generating OSPREY Script..."):
                generated = logic3.create_folders_and_files(
                    json_path, pdb_path,
                    epsilon=epsilon,
                    cpu_cores=int(cores),
                    gpus=int(gpus),
                    streams_per_gpu=int(streams_per_gpu),
                    heap_size=int(heap_size),
                    garbage_size=int(garbage_size),
                    target_residue=target_res,
                    amino_acids=mutations
                )
            st.session_state["osprey_scripts"] = generated

            with st.spinner("Running OSPREY (this may take a while)..."):
                run_results = logic3.run_osprey_scripts(generated)

            # ── CHANGED: save each run's stdout into its residue folder ──
            for r in run_results:
                if r.get("stdout"):
                    out_dir = os.path.dirname(os.path.abspath(r["script"]))
                    with open(os.path.join(out_dir, "osprey_stdout.txt"), "w", encoding="utf-8") as f:
                        f.write(r["stdout"])

            st.session_state["step3_result"] = {
                "target_res":  target_res,
                "generated":   generated,
                "run_results": run_results,
            }
            st.session_state["step3_done"] = all(
                ("error" not in r) and r["returncode"] == 0 for r in run_results
            )
            # Keep a running list of residues that finished successfully so results
            # from every scan can be viewed together in Step 4.
            if st.session_state["step3_done"]:
                done = list(st.session_state.get("scanned_residues", []))
                if target_res not in done:
                    done.append(target_res)
                st.session_state["scanned_residues"] = done
                st.rerun()
                 
        except ValueError as e:
            st.error(f"❌ Setup error: {e}")
            st.stop()

    # ── Show the latest run (persists when you leave and return) ──────────
    res3 = st.session_state.get("step3_result")
    if res3:
        st.success(f"{len(res3['generated'])} OSPREY script(s) generated!")
        for f in res3["generated"]:
            st.write(f"• {f}")

        for r in res3["run_results"]:
            if "error" in r:
                st.error(f"❌ {r['script']}: {r['error']}")
            elif r["returncode"] != 0:
                st.error(f"❌ {r['script']} failed:\n{r['stderr']}")
            else:
                st.success(f"✅ {r['script']} completed")
                if r["stdout"]:
                    with st.expander("View output"):
                        st.code(r["stdout"])

    if st.session_state.get("step3_done", False):
        st.info(f"✅ {target_res} finished. Choose another residue above to scan it too, "
                "or go to Results to see everything scanned so far.")
    next_step_button("4. Results & Analytics", "Click to head to Results & Analytics",
                     ready=bool(st.session_state.get("scanned_residues")))

# ═══════════════════════════════════════════════════════════════════════════════
# Step 4: Results & Analytics
# ═══════════════════════════════════════════════════════════════════════════════
elif workflow_step == "4. Results & Analytics":
    st.header("Step 4: Redesign Results")

    osprey_scripts = st.session_state.get("osprey_scripts", [])
    results_shown = False

    step3_data = st.session_state.get("step3_result")

    if not step3_data or "run_results" not in step3_data:
        st.warning("No OSPREY simulation results found. Please run Step 3 first.")
    else:
        # Pass stdout execution results directly to logic4
        base_dir = st.session_state.get("osprey_dir")
        # ── CHANGED: only include residues scanned in this session ──
        pivot_delta, pivot_annot, summary_table, error = logic4.analyze_results(
            base_dir, residues=st.session_state.get("scanned_residues")
        )

    
        if error:
            st.error(f"Analysis Error: {error}")
        else:
            import matplotlib
            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors
            import numpy as np

            # --- Top metrics row ---
            st.subheader("📊 Affinity Summary")
            col1, col2, col3 = st.columns(3)
            try:
                best_val = pivot_delta.max().max()
                best_mut = pivot_delta.max().idxmax()
                worst_val = pivot_delta.min().min()
                total_runs = pivot_delta.count().sum()
                col1.metric("Top ΔK* Improvement", f"+{best_val:.3f}", f"Mutation: {best_mut}")
                col2.metric("Largest Decrease", f"{worst_val:.3f}")
                col3.metric("Total Sequences Scored", int(total_runs))
            except Exception:
                st.info("Insufficient data for metrics.")

            st.divider()

            # --- Panel A: Heatmap ---
            st.subheader("Panel A — CaM Residue Mutations (ΔK* vs Wild-Type)")
            st.caption("Each cell shows the single-letter amino acid tested. Color = ΔK* (purple = improved affinity, white = neutral).")

            n_rows, n_cols = pivot_delta.shape

            fig_h = max(4, n_rows * 0.55)
            fig_w = max(5, n_cols * 1.1)

            fig, ax = plt.subplots(figsize=(fig_w, fig_h))

            # Color map matching the manuscript (white → purple)
            cmap = matplotlib.colormaps.get_cmap('Purples')

            # Mask NaN cells
            data_vals = pivot_delta.values.astype(float)
            annot_vals = pivot_annot.values

            # Normalize color scale
            valid = data_vals[~np.isnan(data_vals)]
            vmin = float(np.min(valid)) if len(valid) else -0.01
            vmax = float(np.max(valid)) if len(valid) else 0.01
            # TwoSlopeNorm requires vmin < vcenter < vmax strictly
            if vmin >= 0:
                vmin = -0.01
            if vmax <= 0:
                vmax = 0.01
            norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

            # Draw cells manually for full control
            for r in range(n_rows):
                for c in range(n_cols):
                    val = data_vals[r, c]
                    aa  = annot_vals[r, c]
                    if np.isnan(val) or (isinstance(aa, float) and np.isnan(aa)):
                        facecolor = '#f0f0f0'
                        text = ''
                    else:
                        facecolor = cmap(norm(val))
                        text = str(aa) if aa else ''

                    rect = plt.Rectangle([c, r], 1, 1, facecolor=facecolor, edgecolor='white', linewidth=1.5)
                    ax.add_patch(rect)
                    if text:
                        brightness = 0.299*facecolor[0] + 0.587*facecolor[1] + 0.114*facecolor[2]
                        txt_color = 'white' if brightness < 0.55 else '#222222'
                        ax.text(c + 0.5, r + 0.5, text, ha='center', va='center',
                                fontsize=10, fontweight='bold', color=txt_color, fontfamily='monospace')

            # Axes formatting
            ax.set_xlim(0, n_cols)
            ax.set_ylim(0, n_rows)
            ax.set_xticks([c + 0.5 for c in range(n_cols)])
            ax.set_xticklabels(pivot_delta.columns, rotation=35, ha='right', fontsize=9)
            ax.set_yticks([r + 0.5 for r in range(n_rows)])
            ax.set_yticklabels(pivot_delta.index, fontsize=9)
            ax.set_xlabel("Model / Run", fontsize=10, labelpad=8)
            ax.set_ylabel("WT Residue", fontsize=10, labelpad=8)
            ax.tick_params(length=0)
            for spine in ax.spines.values():
                spine.set_visible(False)

            # Colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
            cbar.set_label("ΔK* (log10)", fontsize=9)
            cbar.ax.tick_params(labelsize=8)

            fig.tight_layout()
            st.pyplot(fig)
            plt.close(fig)

            st.divider()

            # --- Summary table & missing AA inspection ---
            st.subheader("🔬 Detailed Sequence Scores")

            AMINO_ACIDS = [
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
            ]

            # Try parsing for missing amino acid tab if standard format exists
            has_parsed_missing = False
            try:
                raw = summary_table.copy()
                raw[["Chain", "Position", "WT"]] = raw["Residue Site"].str.extract(r"([A-Za-z])(\d+)=(\w+)")
                raw["Position"] = pd.to_numeric(raw["Position"], errors="coerce")

                tidy_summary = (
                    raw.assign(pair=raw["Mutations & Log10 Scores"].str.findall(r"(\w+)\s*\(([\d.\-eE]+)\)"))
                    .explode("pair")
                    .dropna(subset=["pair"])
                )
                tidy_summary[["Mutation", "Log10 Score"]] = pd.DataFrame(
                    tidy_summary["pair"].tolist(), index=tidy_summary.index
                )
                tidy_summary["Log10 Score"] = tidy_summary["Log10 Score"].astype(float)
                tidy_summary["Mutation"] = tidy_summary["Mutation"].str.upper()

                tidy_summary = (
                    tidy_summary
                    .drop(columns=["pair", "Mutations & Log10 Scores"])
                    .drop_duplicates(subset=["Residue Site", "Mutation"])
                    .sort_values(["Chain", "Position", "Log10 Score"], ascending=[True, True, False])
                    .reset_index(drop=True)
                )
                tidy_summary = tidy_summary[["Residue Site", "Chain", "Position", "WT", "Mutation", "Log10 Score"]]

                site_info = tidy_summary.drop_duplicates("Residue Site")[["Residue Site", "Chain", "Position", "WT"]]
                tested_by_site = tidy_summary.groupby("Residue Site")["Mutation"].apply(set)

                missing_by_site = site_info.copy()
                missing_by_site["Missing Amino Acids"] = missing_by_site["Residue Site"].map(
                    lambda s: sorted(set(AMINO_ACIDS) - tested_by_site.get(s, set()))
                )
                missing_by_site["# Missing"] = missing_by_site["Missing Amino Acids"].str.len()
                missing_by_site["Missing Amino Acids"] = missing_by_site["Missing Amino Acids"].apply(", ".join)
                missing_by_site = missing_by_site.sort_values(["Chain", "Position"]).reset_index(drop=True)

                missing_globally = sorted(set(AMINO_ACIDS) - set(tidy_summary["Mutation"]))
                has_parsed_missing = True
            except Exception:
                has_parsed_missing = False

            if has_parsed_missing:
                tab1, tab2, tab3 = st.tabs(["Summary Table", "Raw ΔK* Matrix", "Missing Amino Acids"])
                with tab1:
                    st.write("All tested mutations and their absolute K* scores:")
                    st.dataframe(
                        tidy_summary,
                        use_container_width=True,
                        hide_index=True,
                        column_config={"Log10 Score": st.column_config.NumberColumn(format="%.3f")},
                    )
                with tab2:
                    st.write("ΔK* matrix (0.000 = Wild-Type reference):")
                    st.dataframe(pivot_delta.style.format("{:.3f}").background_gradient(cmap='Purples'), use_container_width=True)
                with tab3:
                    st.write("Standard amino acids not represented in the summary table:")
                    if missing_globally:
                        st.warning(f"Never tested at any site: {', '.join(missing_globally)}")
                    else:
                        st.success("All 20 amino acids appear at least once across the dataset.")

                    st.write("Missing amino acids per residue site:")
                    st.dataframe(missing_by_site, use_container_width=True, hide_index=True)
            else:
                tab1, tab2 = st.tabs(["Summary Table", "Raw ΔK* Matrix"])
                with tab1:
                    st.write("All tested mutations and their absolute K* scores:")
                    st.dataframe(summary_table, use_container_width=True, hide_index=True)
                with tab2:
                    st.write("ΔK* matrix (0.000 = Wild-Type reference):")
                    st.dataframe(pivot_delta.style.format("{:.3f}").background_gradient(cmap='Purples'), use_container_width=True)

            # --- Export ---
            st.divider()
            st.subheader("📥 Export Data")
            csv = pivot_delta.to_csv().encode('utf-8')
            st.download_button(
                label="Download Results as CSV",
                data=csv,
                file_name="osprey_redesign_results.csv",
                mime="text/csv",
            )
            results_shown = True
    next_step_button("3. OSPREY Execution", "Scan another residue", ready=results_shown)
# ── Sidebar footer ────────────────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.write("**Software Team (Davis Lab)**")