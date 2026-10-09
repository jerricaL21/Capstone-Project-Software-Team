import streamlit as st
import pandas as pd
import time
import os
import importlib
from install_pymol import ensure_pymol

# ── Ensure PyMOL is available in this interpreter (cached so it only
#    checks/installs once per session, not on every Streamlit rerun) ──
@st.cache_resource
def _pymol_ready():
    return ensure_pymol()

pymol_ok = _pymol_ready()
if not pymol_ok:
    st.sidebar.error("⚠️ PyMOL not available — Step 5 will not work.")

# Importing the script files
import logic
import logic2
import logic3
import logic4
import logic5
import logic6
import logic7

importlib.reload(logic)
importlib.reload(logic2)
importlib.reload(logic3)
importlib.reload(logic4)
importlib.reload(logic5)
importlib.reload(logic6)
importlib.reload(logic7)

# Used in Step 3 to figure out which three-letter amino acid code corresponds to wild-type residue code
AA_1TO3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}

#Page Configuration
st.set_page_config(page_title="Davis Lab Redesign", layout="wide")

#UI Header
st.title("CaM-RyR2 Protein Redesign Portal")
st.markdown("""
**Status:** Unified Workflow v1.0 (Draft)  
*Meeting ISO/IEC 25010:2023 Efficiency Standards*
""")

st.sidebar.header("Pipeline Controls")
workflow_step = st.sidebar.radio("Navigate Workflow", 
    ["1. Data Input", "2. JSON Processing", "3. OSPREY Execution", "4. Results & Analytics", "5. PyMOL Redesign", "6. Molecular Dynamics"],
    key="workflow_step")

# ── "Next step" button shown at the bottom of each page ──────────────────────
def _go_to_step(target):
    st.session_state["workflow_step"] = target

def next_step_button(target, label, ready = True):
    if not ready:
        return
    st.markdown("---")
    st.button(f"{label}  ➡️", key=f"next_to_{target}", on_click=_go_to_step,
              args=(target,), type="primary", use_container_width=True)

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Data Input
# ═══════════════════════════════════════════════════════════════════════════════
if workflow_step == "1. Data Input":
    st.header("Step 1: Structural Data Input")

    uploaded_pdb = st.file_uploader("Upload Target PDB File", type=['pdb'], key="pdb_uploader")

    # ── A different PDB was uploaded: reset everything, then load the new file ──
    if uploaded_pdb is not None and uploaded_pdb.name != st.session_state.get("pdb_name"):
        _PROTECTED = {"workflow_step", "pdb_uploader"}
        for k in list(st.session_state.keys()):
            if k not in _PROTECTED:
                del st.session_state[k]

        os.makedirs("uploaded_pdbs", exist_ok=True)
        new_path = os.path.abspath(os.path.join("uploaded_pdbs", uploaded_pdb.name))
        with open(new_path, "wb") as f:
            f.write(uploaded_pdb.getvalue())

        st.session_state["pdb_name"]  = uploaded_pdb.name
        st.session_state["orig_path"] = new_path   # untouched upload; Run Analysis always starts from this
        st.session_state["path"]      = new_path   # later steps read this (updated after repair/protonation)
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
                    if info["is_calmodulin"]:
                        if info.get("identity") is not None:
                            identity_msg = f" ({info['identity']}% identity to human CaM, {info['length']} residues)"
                        else:
                            identity_msg = f" ({info['length']} residues)"
                        st.success(f"Chain **{chain_id}** → CALMODULIN{identity_msg}")
                    else:
                        st.info(f"Chain **{chain_id}** → Ligand / peptide ({info['length']} residues)")

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
                    "Select which chain is the **calmodulin** and which is the "
                    "**peptide/ligand (residue) chain**."
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
                        "Calmodulin chain",
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
                        "Peptide / ligand chain",
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
                st.session_state["json_path"]       = json_path
                st.session_state["step1_done_for"]  = pdb_name
                st.session_state["analysis_result"] = {"repair": repair_info, "json_path": json_path}
                st.session_state.pop("step2_done", None)
                st.session_state.pop("step3_done", None)
                st.session_state.pop("step3_result", None)

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
# Step 2: JSON Processing  (unchanged)
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
# Step 3: OSPREY Execution  (unchanged)
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
        target_res = st.selectbox(
            "Select Residue to Mutate", residues,
            index=residues.index(prev_res) if prev_res in residues else 0,
        )

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

            st.session_state["step3_result"] = {
                "target_res":  target_res,
                "generated":   generated,
                "run_results": run_results,
            }
            st.session_state["step3_done"] = all(
                ("error" not in r) and r["returncode"] == 0 for r in run_results
            )    
                 
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

    next_step_button("4. Results & Analytics", "Click to head to Results & Analytics", ready=st.session_state.get("step3_done", False))

# ═══════════════════════════════════════════════════════════════════════════════
# Step 4: Results & Analytics  (unchanged)
# ═══════════════════════════════════════════════════════════════════════════════
#Step 4: Results
elif workflow_step == "4. Results & Analytics":
    st.header("Step 4: Redesign Results")

    osprey_scripts = st.session_state.get("osprey_scripts", [])
    results_shown = False

    if not osprey_scripts:
        st.warning("No OSPREY scripts found. Please complete Step 3 first.")
    else:
        base_dir = os.path.dirname(os.path.dirname(osprey_scripts[0]))

        pivot_delta, pivot_annot, summary_table, error = logic4.analyze_results(base_dir)

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

            # --- Summary table ---
            st.subheader("🔬 Detailed Sequence Scores")
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

    next_step_button("5. PyMOL Redesign", "Click to head to PyMOL Redesign", ready=results_shown)
# ═══════════════════════════════════════════════════════════════════════════════
# Step 5: PyMOL Redesign
# ═══════════════════════════════════════════════════════════════════════════════
elif workflow_step == "5. PyMOL Redesign":
    st.header("Step 5: Multi-Site Redesign (PyMOL)")
    st.markdown(
        "This step combines the best-scoring mutation from **selected** residue "
        "sites you've scanned in Step 3 into a single structure."
    )

    osprey_scripts = st.session_state.get("osprey_scripts", [])
    pdb_for_pymol  = st.session_state.get("path", None)

    if not osprey_scripts or not pdb_for_pymol:
        st.warning("No OSPREY scans found yet. Please complete Step 3 first.")
        st.stop()

    base_dir = os.path.dirname(os.path.dirname(osprey_scripts[0]))

    if st.button("🔄 Check Scanned Sites"):
        st.session_state.pop("pymol_pivot_delta", None)

    if "pymol_pivot_delta" not in st.session_state:
        with st.spinner("Aggregating mutation-scan results..."):
            pivot_delta, _pivot_annot, _summary, error = logic4.analyze_results(base_dir)
        if error:
            st.error(f"Analysis Error: {error}")
            st.stop()
        st.session_state["pymol_pivot_delta"] = pivot_delta
    else:
        pivot_delta = st.session_state["pymol_pivot_delta"]

    if pivot_delta is None or pivot_delta.empty:
        st.warning("No residue sites found in analysis results. Run Step 3 first.")
        st.stop()

    available_sites = list(pivot_delta.index)

    # ── Allow user to pick specific residue sites via multiselect ─────────────
    st.subheader("🎯 Select Residues to Include")
    selected_sites = st.multiselect(
        "Choose specific residue sites to send to PyMOL:",
        options=available_sites,
        default=available_sites,  # Defaults to selecting all available sites
        help="Select one or more residue sites to mutate."
    )

    if not selected_sites:
        st.info("Please select at least one residue site above to proceed.")
        st.stop()

    # ── Preview the best mutation per selected site ───────────────────────────
    try:
        preview_mutations = logic6.select_best_mutations(pivot_delta, selected_sites=selected_sites)
    except ValueError as e:
        st.error(f"❌ {e}")
        st.stop()

    st.subheader(f"🏆 Best Mutation for {len(preview_mutations)} Selected Site(s)")
    st.dataframe(
        pd.DataFrame(preview_mutations)[["site", "chain", "resnum", "wt_aa3", "mutant_aa3", "delta"]]
          .rename(columns={
              "site": "Site", "chain": "Chain", "resnum": "Residue #",
              "wt_aa3": "Wild-Type", "mutant_aa3": "Best Mutation", "delta": "ΔK*",
          }),
        use_container_width=True, hide_index=True,
    )

    st.caption(f"Structure to be mutated: `{pdb_for_pymol}`")

    if st.button("🧬 Generate Redesigned Structure in PyMOL"):
        try:
            with st.spinner("Running PyMOL mutagenesis..."):
                output_path, mutations, applied, errors = logic6.generate_redesigned_structure(
                    pdb_for_pymol, pivot_delta, selected_sites=selected_sites
                )
            st.session_state["pymol_output_path"] = output_path

            if applied:
                st.success(f"✅ Applied {len(applied)} mutation(s):")
                for a in applied:
                    st.write(f"• {a}")
            if errors:
                st.warning("⚠️ Some sites could not be mutated:")
                for err in errors:
                    st.write(f"• {err}")

            st.success(f"Redesigned PDB saved to: `{output_path}`")

        except ImportError as e:
            st.error(f"❌ {e}")
        except ValueError as e:
            st.error(f"❌ {e}")

    output_path = st.session_state.get("pymol_output_path")
    if output_path and os.path.exists(output_path):
        with open(output_path, "rb") as f:
            st.download_button(
                label="📥 Download Redesigned PDB",
                data=f.read(),
                file_name=os.path.basename(output_path),
                mime="chemical/x-pdb",
            )
    next_step_button("6. Molecular Dynamics", "Click to head to Molecular Dynamics (NAMD)",
                     ready=bool(output_path and os.path.exists(output_path)))
# ═══════════════════════════════════════════════════════════════════════════════
# Step 6: Molecular Dynamics
# ═══════════════════════════════════════════════════════════════════════════════
elif workflow_step == "6. Molecular Dynamics":
    st.header("Step 6: Molecular Dynamics — NAMD Package Generator")

    st.markdown("""\
    This step generates a **ready-to-run, explicit-solvent NAMD package** for your PDB structure.
    The portal does not execute NAMD directly (simulations take hours to days on a GPU), but packages
    everything your group needs — a solvated CHARMM PSF/PDB, a staged 4-part NAMD protocol
    (minimize → heat → equilibrate → production), a launcher script, and MM/PBSA scaffolding —
    so you can run it on your own GPU workstation/cluster **without VMD**.

    > Fig. 2B tracks the N-domain ↔ C-domain distance over time.  
    > - **wtCaM** → stays in *annealed* state (~0.5–1 nm) → normal regulation  
    > - **RCaM1** → unlocks (~1.5–3 nm), bends RyR2 peptide → Ca²⁺ leak ↑  
    > - **RCaM2** → stays locked and annealed → Ca²⁺ leak ↓ (therapeutic goal)
    """)

    st.divider()

    # A: PDB source (PATH FIX: Checks Step 5 output first, falls back to Step 1)
    st.subheader("A. Structure Source")
    pdb_from_session = st.session_state.get("pymol_output_path") or st.session_state.get("path", None)
    md_pdb_path = None

    if pdb_from_session and os.path.exists(pdb_from_session):
        source_label = "Step 5 (PyMOL Redesign)" if pdb_from_session == st.session_state.get("pymol_output_path") else "Step 1"
        use_session = st.checkbox(
            f"Use PDB from {source_label}: `{os.path.basename(pdb_from_session)}`", value=True)
        if use_session:
            md_pdb_path = pdb_from_session
    else:
        use_session = False

    if not use_session or md_pdb_path is None:
        uploaded_md_pdb = st.file_uploader("Upload PDB file for MD", type=["pdb"], key="md_pdb")
        if uploaded_md_pdb:
            os.makedirs("uploaded_pdbs", exist_ok=True)
            md_pdb_path = os.path.abspath(os.path.join("uploaded_pdbs", uploaded_md_pdb.name))
            with open(md_pdb_path, "wb") as fh:
                fh.write(uploaded_md_pdb.getvalue())
            st.success(f"Uploaded: `{md_pdb_path}`")
    st.divider()

    # B: Parameters
    st.subheader("B. Simulation Parameters")
    st.caption("Explicit TIP3P solvent only — matches the manuscript's GROMACS methodology. "
               "(A previous version of this tool offered implicit GB solvent; that option was "
               "removed because it does not reproduce the paper's setup.)")
    col1, col2 = st.columns(2)
    with col1:
        temperature  = st.number_input("Temperature (K)", 270.0, 400.0, 297.15, step=5.0,
                                        help="297.15 K = 24 °C — matches the SI Appendix's MD Methods "
                                             "(not physiological 310 K).")
        production_ns = st.number_input("PRODUCTION stage length (ns)", 0.1, 500.0, 100.0, step=0.1,
                                        help="SI Appendix: production run of at least 100 ns. This is "
                                             "on top of the fixed minimize/heat/equilibration stages, "
                                             "not instead of them.")
        timestep_fs  = float(st.selectbox("Timestep (fs)", [1.0, 2.0], index=1))
    with col2:
        nonbonded_cutoff = st.number_input("Non-bonded cutoff (Å)", 8.0, 20.0, 12.0, step=1.0)
        gpu_enabled      = st.checkbox("Enable CUDA GPU acceleration", value=True,
                                        help="Requires NAMD CUDA build — recommended for Windows + RTX/GTX")
        padding_nm       = st.number_input("Water box padding (nm)", 0.8, 3.0, 1.0, step=0.1,
                                        help="SI Appendix used a minimum 1.0 nm (10 Å) protein-to-box-edge "
                                             "distance. Increase if your variant undergoes large conformational "
                                             "excursions (e.g. RCaM1-like unlocking) to avoid self-interaction "
                                             "across periodic images.")
        ionic_strength   = st.number_input("Ionic strength (M KCl)", 0.0, 1.0, 0.15, step=0.05,
                                        help="SI Appendix used KCl, not NaCl, at 0.15 M.")
        protein_chain    = st.text_input("Protein chain ID", value="A")
        ligand_chain     = st.text_input("Peptide/ligand chain ID", value="B")

    st.session_state["md_timestep_fs"] = timestep_fs

    st.divider()

    # C: Generate (PATH & EXECUTION ROUTING FIX: Routes directly via logic7)
    st.subheader("C. Generate NAMD Package")

    if md_pdb_path is None:
        st.info("Upload or select a PDB file in section A to enable package generation.")
    else:
        st.write(f"Ready to package: **`{os.path.basename(md_pdb_path)}`**")
        st.caption("This step solvates the structure (TIP3P + ions), builds a CHARMM PSF, and writes a "
                   "4-stage NAMD protocol (minimize → heat → restrained NPT equilibration → unrestrained "
                   "NPT production) plus MM/PBSA scaffolding. Solvation runs here; NAMD itself does not — "
                   "run the generated package on your own GPU workstation/cluster.")
        if st.button("Generate NAMD Package"):
            with st.spinner("Solvating system, building PSF, and writing staged NAMD configs... this can take a minute or two."):
                out_dir = os.path.dirname(os.path.abspath(md_pdb_path))
                try:
                    # Routing execution to teammates' logic7
                    md_func = getattr(logic7, "run_md_pipeline_from_pymol", None) or getattr(logic7, "build_full_package", None)
                    if md_func is None:
                        raise AttributeError("Neither 'run_md_pipeline_from_pymol' nor 'build_full_package' was found in logic7.py")

                    result = md_func(
                        pdb_path=md_pdb_path, output_dir=out_dir,
                        padding_nm=padding_nm, ionic_strength_m=ionic_strength,
                        protein_chain=protein_chain, ligand_chain=ligand_chain,
                        temperature=temperature, production_ns=production_ns,
                        timestep_fs=timestep_fs, nonbonded_cutoff=nonbonded_cutoff,
                        gpu_enabled=gpu_enabled,
                    )
                except Exception as e:
                    st.error(f"MD packaging failed: {e}")
                    st.stop()
            st.success(result["summary"])
            if result.get("manifest", {}).get("warnings"):
                for w in result["manifest"]["warnings"]:
                    st.warning(w)
            st.session_state["namd_result"] = result
            st.session_state["md_dcd_freq"] = result.get("namd_package", {}).get("production_dcd_freq", 5000)

            namd_dir = result["namd_package"]["namd_dir"]
            labels = {"README": result["namd_package"]["readme_path"]}
            for stage in result["namd_package"]["stages"]:
                labels[f"NAMD config: {stage}"] = os.path.join(namd_dir, stage)
            labels["Analysis script (analyze_distance.py)"] = result["namd_package"]["analysis_path"]
            labels["Launcher (run_all_stages.sh)"] = result["namd_package"]["launcher_path"]
            labels["MM/PBSA README"] = os.path.join(result["mmpbsa"]["mmpbsa_dir"], "README.txt")

            for label, fpath in labels.items():
                if os.path.exists(fpath):
                    with st.expander(f"View: {label}"):
                        with open(fpath) as fh:
                            st.code(fh.read(), language="bash")

            st.divider()
            if "zip_path" in result and os.path.exists(result["zip_path"]):
                with open(result["zip_path"], "rb") as zf:
                    st.download_button(
                        label="Download complete NAMD package (.zip)",
                        data=zf.read(),
                        file_name=os.path.basename(result["zip_path"]),
                        mime="application/zip",
                    )

    st.divider()

    # D: Upload & plot MD results
    st.subheader("D. Visualise MD Results — N/C Domain Distance")
    st.markdown("""\
    After running NAMD and the `analyze_distance.py` script, upload `domain_distance.dat`
    to plot the N-domain / C-domain distance over time (replicates manuscript Fig. 2B).
    """)

    dat_file = st.file_uploader("Upload domain_distance.dat", type=["dat","txt","csv"], key="dat_upload")
    if dat_file:
        import matplotlib.pyplot as plt

        dat_text = dat_file.read().decode("utf-8")
        rows     = logic6.parse_domain_distance_dat(dat_text)

        if not rows:
            st.error("Could not parse the file. Expected two columns: frame and distance (nm).")
        else:
            df = pd.DataFrame(rows)
            ts = st.session_state.get("md_timestep_fs", 2.0)
            dcd_freq = st.session_state.get("md_dcd_freq", 5000)  # falls back to this session's default if package wasn't generated in this run
            df["time_ns"] = df["frame"] * dcd_freq * ts / 1_000_000

            traj_label = st.text_input("Trajectory label (e.g. wtCaM, RCaM1, RCaM2)", value="My CaM")

            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(df["time_ns"], df["distance_nm"], linewidth=1.0,
                    label=traj_label, color="#6a0dad")
            ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8,
                       label="Annealed threshold (~1 nm)")
            ax.set_xlabel("Time (ns)", fontsize=11)
            ax.set_ylabel("N-C Domain Distance (nm)", fontsize=11)
            ax.set_title("CaM N-domain / C-domain Distance Over Time", fontsize=12)
            ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
            fig.tight_layout(); st.pyplot(fig); plt.close(fig)

            c1, c2, c3 = st.columns(3)
            c1.metric("Mean distance", f"{df['distance_nm'].mean():.2f} nm")
            c2.metric("Min distance",  f"{df['distance_nm'].min():.2f} nm")
            c3.metric("Max distance",  f"{df['distance_nm'].max():.2f} nm")

            annealed_pct = (df["distance_nm"] < 1.0).mean() * 100
            st.caption(
                f"Fraction of frames in annealed state (<1 nm): **{annealed_pct:.1f}%**  "
                "(wtCaM/RCaM2 ≈ high, RCaM1 ≈ low)"
            )

    st.divider()

    with st.expander("Quick-Start Reference — NAMD on Windows (no VMD needed)"):
        st.markdown("""\
**Required software (all free)**
- [NAMD 2.14 multicore](https://www.ks.uiuc.edu/Research/namd/) — Win64, no GPU version works fine
- [CHARMM36 force field](http://mackerell.umaryland.edu/charmm_ff.shtml) — `toppar_c36_jul22.tgz`
- parmed + MDAnalysis — `pip install parmed MDAnalysis`

**Steps after downloading the .zip**
1. Extract the zip. Copy your PDB into the folder. Copy `top_all36_prot.rtf` and `par_all36_prot.prm` from the CHARMM36 download into the same folder.
2. Generate the topology (replaces VMD psfgen):
   `python generate_psf.py`
3. Edit the `.namd` file — change `coordinates` to point to `{pdb_name}_prep.pdb`
4. Launch NAMD:
   `namd2.exe +p4 {pdb_name}.namd > {pdb_name}_md.log`
5. After simulation, compute the domain distance:
   `python analyze_distance.py`
6. Upload `domain_distance.dat` to section D above.

**Interpreting the distance plot**

| CaM variant | Distance behavior | Biological meaning |
|-------------|-------------------|-------------------|
| wtCaM | Stays low (<1 nm), "annealed" | Normal RyR2 regulation |
| RCaM1 | Rises high (>1.5 nm), "unlocked" | Bends RyR2 peptide → Ca²⁺ leak ↑ |
| RCaM2 | Stays low, "locked annealed" | High affinity + straight peptide → Ca²⁺ leak ↓ |
        """)


# ── Sidebar footer ────────────────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.write("**Software Team (Davis Lab)**")