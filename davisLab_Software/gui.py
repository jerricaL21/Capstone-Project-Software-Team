import streamlit as st
import pandas as pd
import time
import os
import importlib

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
    ["1. Data Input", "2. JSON Processing", "3. OSPREY Execution", "4. Results & Analytics", "5. Mutant Prediction"])

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Data Input
# ═══════════════════════════════════════════════════════════════════════════════
if workflow_step == "1. Data Input":
    st.header("Step 1: Structural Data Input")
    
    uploaded_pdb = st.file_uploader("Upload Target PDB File", type=['pdb'])
    
    if uploaded_pdb is not None:
        
        os.makedirs("uploaded_pdbs", exist_ok=True)
        path = os.path.abspath(os.path.join("uploaded_pdbs", uploaded_pdb.name))
        
        with open(path, "wb") as f:
            f.write(uploaded_pdb.getvalue())
        
        st.session_state["path"] = path

        pdb_bytes = uploaded_pdb.getvalue().decode("utf-8")
        
        st.success(f"Successfully loaded: {uploaded_pdb.name}")
        st.write(f"Saved file path: {path}")
        
        with st.expander("📄 View PDB File Header"):
            st.code(pdb_bytes[:500] + "...", language="text")
        
        st.warning("Note: To view the 3D structure, you will need to install 'stmol'.")

        # ── Chain identification — only runs once per uploaded file ───────
        already_identified = (
            "chain_results" in st.session_state and
            st.session_state.get("identified_for") == uploaded_pdb.name
        )

        if not already_identified:
            with st.spinner("Identifying chains..."):
                try:
                    chain_results = logic.identify_calmodulin_chain(path)
                    st.session_state["chain_results"]    = chain_results
                    st.session_state["identified_for"]   = uploaded_pdb.name
                except Exception as e:
                    st.error(f"Chain identification failed: {e}")
                    chain_results = {}
        else:
            chain_results = st.session_state["chain_results"]

        # ── Display chain identification results ──────────────────────────
        if chain_results:
            st.subheader("🔍 Chain Identification")

            cam_chain    = None
            ligand_chain = None

            for chain_id, info in chain_results.items():
                if info["is_calmodulin"]:
                    cam_chain = chain_id
                    st.success(
                        f"Chain **{chain_id}** → CALMODULIN "
                        f"({info['identity']}% identity to human CaM, {info['length']} residues)"
                    )
                else:
                    ligand_chain = chain_id
                    st.info(
                        f"Chain **{chain_id}** → Ligand / peptide "
                        f"({info['identity']}% identity to human CaM, {info['length']} residues)"
                    )

            if not cam_chain:
                st.warning(
                    "⚠️ No chain matched calmodulin (≥80% identity). "
                    "Proceed to Step 2 and enter residue ranges manually."
                )

            st.session_state["cam_chain"]    = cam_chain
            st.session_state["ligand_chain"] = ligand_chain

        # ── PDB Inspection (NEW) ──────────────────────────────────────────
        # Runs logic5.inspect_pdb() and shows the user what problems were
        # found before the repair step.
        already_inspected = (
            "inspection_report" in st.session_state and
            st.session_state.get("inspected_for") == uploaded_pdb.name
        )

        if not already_inspected:
            with st.spinner("Inspecting PDB for structural completeness..."):
                try:
                    report = logic5.get_inspection_report(path)
                    st.session_state["inspection_report"] = report
                    st.session_state["inspected_for"]     = uploaded_pdb.name
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
            if report["has_missing_residues"]:
                st.warning(
                    "⚠️ **Missing residues detected.** Automatic remodelling is currently "
                    "disabled — the pipeline will proceed with the original coordinates. "
                    "Side chains will still be repaired if needed. "
                    "Boltz-2 residue remodelling can be enabled in `logic5.py` when ready."
                )
            elif report["has_missing_sidechains"]:
                st.info(
                    "🔧 **Repair strategy: PDBFixer** — Only side-chain atoms "
                    "are missing. Fast in-process repair, no GPU needed."
                )
            else:
                st.success("✅ PDB is complete — no repair needed before analysis.")

        # ── MSA server option (only relevant if Boltz-2 will be used) ────
        use_msa = True
        if report and report["has_missing_residues"]:
            use_msa = st.checkbox(
                "Use ColabFold MSA server (requires internet connection)",
                value=True,
                help="Boltz-2 uses multiple sequence alignments to improve "
                     "prediction quality. Uncheck if you are air-gapped or "
                     "have pre-generated MSAs."
            )
            st.session_state["use_msa_server"] = use_msa

        # ── Run Analysis button ───────────────────────────────────────────
        if st.button("▶ Run Analysis"):

            cam_chain    = st.session_state.get("cam_chain")
            ligand_chain = st.session_state.get("ligand_chain")

            # Step A: PDB repair (NEW — runs logic5 before prppi)
            if report and (report["has_missing_residues"] or report["has_missing_sidechains"]):
                if not cam_chain or not ligand_chain:
                    st.error(
                        "❌ Chain identification is required before repair. "
                        "Please ensure calmodulin and peptide chains were identified above."
                    )
                    st.stop()

                with st.spinner("Repairing PDB structure..."):
                    try:
                        repaired_path, strategy = logic5.prepare_pdb(
                            pdb_path         = path,
                            cam_chain_id     = cam_chain,
                            peptide_chain_id = ligand_chain,
                            use_msa_server   = st.session_state.get("use_msa_server", True),
                        )
                        st.session_state["path"] = repaired_path

                        strategy_labels = {
                            "none_needed":              "No repair was needed.",
                            "sidechain_only":           "Side chains rebuilt with PDBFixer.",
                            "boltz_repair":             "Missing residues reconstructed with Boltz-2.",
                            "missing_residues_skipped": "Missing residues detected but remodelling is currently disabled — proceeding with original coordinates.",
                        }
                        st.success(
                            f"✅ Repair complete: {strategy_labels.get(strategy, strategy)}\n\n"
                            f"Repaired file: `{repaired_path}`"
                        )
                        path = repaired_path

                    except Exception as e:
                        st.error(f"❌ PDB repair failed: {e}")
                        st.stop()
            else:
                st.info("PDB is complete — skipping repair step.")

            # Step B: Fix protonation (existing logic.py step)
            with st.spinner("Fixing PDB protonation..."):
                fixed_pdb = logic.fix_pdb_protonation(path)
                st.session_state["path"] = fixed_pdb

            # Step C: Run prppi (existing logic.py step)
            with st.spinner("Running PRPPI analysis..."):
                json_path = logic.run_prppi(fixed_pdb, cutoff=5.0)

            if json_path:
                st.success(f"Analysis complete! Results saved to: {json_path}")
                st.session_state["json_path"] = json_path

                with st.expander("🔬 View JSON Contents"):
                    with open(json_path) as f:
                        import json as jsonlib
                        json_data = jsonlib.load(f)
                    st.json(json_data)


# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: JSON Processing  (unchanged)
# ═══════════════════════════════════════════════════════════════════════════════
elif workflow_step == "2. JSON Processing":
    st.header("Step 2: JSON Processing")

    json_path = st.session_state.get("json_path", None)

    if json_path:
        st.success(f"JSON ready from Step 1: {json_path}")
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

        with st.form("residue_form"):

            st.markdown("**Calmodulin (Mutant) Range**")
            if cam_options:
                col1, col2 = st.columns(2)
                with col1:
                    mutant_start = st.selectbox(
                        "Mutant Beginning",
                        options=cam_options,
                        index=0
                    )
                with col2:
                    mutant_end = st.selectbox(
                        "Mutant Ending",
                        options=cam_options,
                        index=len(cam_options) - 1
                    )
            else:
                col1, col2 = st.columns(2)
                with col1:
                    mutant_start = st.text_input("Mutant Beginning")
                with col2:
                    mutant_end = st.text_input("Mutant Ending")

            st.markdown("**Ligand Range**")
            if ligand_options:
                col1, col2 = st.columns(2)
                with col1:
                    ligand_start = st.selectbox(
                        "Ligand Beginning",
                        options=ligand_options,
                        index=0
                    )
                with col2:
                    ligand_end = st.selectbox(
                        "Ligand Ending",
                        options=ligand_options,
                        index=len(ligand_options) - 1
                    )
            else:
                col1, col2 = st.columns(2)
                with col1:
                    ligand_start = st.text_input("Ligand Beginning")
                with col2:
                    ligand_end = st.text_input("Ligand Ending")

            submitted = st.form_submit_button("▶ Process JSON(s)")

        if submitted:
            if not all([mutant_start, mutant_end, ligand_start, ligand_end]):
                st.error("❌ Please fill in all four residue range fields before proceeding.")
            else:
                for json_file in json_paths:
                    logic2.add_variables_to_json(json_file, mutant_start, mutant_end, ligand_start, ligand_end)
                st.success("All JSON files processed! Proceed to Step 3.")


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
        target_res = st.selectbox("Select Residue to Mutate", residues)
        all_aa_options = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'GLU', 'TYR', 'ASP', 'ARG', 'ASN', 'CYS', 'GLN', 'GLY', 'HIS', 'LYS', 'PRO', 'SER', 'THR']
        wt_letter = target_res[-1].upper() if target_res else None
        wt_three_letter = AA_1TO3.get(wt_letter)
        if wt_three_letter:
            mutation_options = [aa for aa in all_aa_options if aa != wt_three_letter]
        else:
            mutation_options = all_aa_options
        mutations = st.multiselect("Select Amino Acids", mutation_options)

    with col2:
        st.subheader("Parameters")
        epsilon = st.slider("K* Precision (Epsilon)", 0.01, 1.0, 0.99)
        cores = st.number_input("CPU Cores", 1, 16, 4)
        gpus = st.number_input("GPUs", min_value=0, max_value=8, value=0, step=1)
        if gpus > 0:
            streams_per_gpu = st.number_input("Streams per GPU", min_value=0, max_value=256, value=82, step=1)
        else:
            streams_per_gpu = 0
            st.info("Streams per GPU set to 0 (no GPU detected)")
        heap_size = st.number_input("Heap Size (MiB)", 1000, 500000, 100000)
        garbage_size = st.number_input("Garbage Size (MiB)", 512, 50000, 8192)

    if st.button("Generate & Run OSPREY Script"):

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
            st.success(f"{len(generated)} OSPREY script(s) generated!")
            for f in generated:
                st.write(f"• {f}")
            st.session_state["osprey_scripts"] = generated

            with st.spinner("Running OSPREY (this may take a while)..."):
                run_results = logic3.run_osprey_scripts(generated)

            for r in run_results:
                if "error" in r:
                    st.error(f"❌ {r['script']}: {r['error']}")
                elif r["returncode"] != 0:
                    st.error(f"❌ {r['script']} failed:\n{r['stderr']}")
                else:
                    st.success(f"✅ {r['script']} completed")
                    if r["stdout"]:
                        with st.expander("View output"):
                            st.code(r["stdout"])

        except ValueError as e:
            st.error(f"❌ Setup error: {e}")
            st.stop()


# ═══════════════════════════════════════════════════════════════════════════════
# Step 4: Results & Analytics  (unchanged)
# ═══════════════════════════════════════════════════════════════════════════════
#Step 4: Results
elif workflow_step == "4. Results & Analytics":
    st.header("Step 4: Redesign Results")

    osprey_scripts = st.session_state.get("osprey_scripts", [])

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

# ═══════════════════════════════════════════════════════════════════════════════
# Step 5: Mutant Prediction (NEW)
# Predicts the structure of a user-specified mutant sequence using Boltz-2,
# then feeds it through the full prppi → OSPREY pipeline.
# ═══════════════════════════════════════════════════════════════════════════════
elif workflow_step == "5. Mutant Prediction":
    st.header("Step 5: Mutant Sequence Structure Prediction")
    st.markdown(
        "Predict the 3D structure of a mutant calmodulin–peptide complex using "
        "Boltz-2, then run it through the full pipeline."
    )

    # ── Require a wildtype PDB to use as template ─────────────────────────
    wt_pdb     = st.session_state.get("path")
    cam_chain  = st.session_state.get("cam_chain")
    pep_chain  = st.session_state.get("ligand_chain")

    if not wt_pdb or not os.path.exists(wt_pdb):
        st.warning("⚠️ Please complete Step 1 first to provide a wildtype structure.")
        st.stop()

    if not cam_chain or not pep_chain:
        st.warning("⚠️ Chain identification not found. Please re-run Step 1.")
        st.stop()

    st.info(
        f"Wildtype structure: `{os.path.basename(wt_pdb)}`  \n"
        f"CaM chain: **{cam_chain}** | Peptide chain: **{pep_chain}**"
    )

    # ── Extract wildtype sequences for display ────────────────────────────
    try:
        wt_cam_seq = logic5.extract_sequence_from_pdb(wt_pdb, cam_chain)
        wt_pep_seq = logic5.extract_sequence_from_pdb(wt_pdb, pep_chain)
    except Exception as e:
        st.error(f"Could not read sequences from wildtype PDB: {e}")
        st.stop()

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Calmodulin (CaM) Sequence**")
        with st.expander("View wildtype CaM sequence"):
            st.code(wt_cam_seq, language="text")
        mutant_cam = st.text_area(
            "Mutant CaM sequence (leave blank to keep wildtype)",
            height=100,
            placeholder=f"Paste mutant sequence here, or leave blank...",
        ).strip()

    with col2:
        st.markdown("**Peptide Sequence**")
        with st.expander("View wildtype peptide sequence"):
            st.code(wt_pep_seq, language="text")
        mutant_pep = st.text_area(
            "Mutant peptide sequence (leave blank to keep wildtype)",
            height=100,
            placeholder=f"Paste mutant sequence here, or leave blank...",
        ).strip()

    use_msa = st.checkbox(
        "Use ColabFold MSA server (requires internet)",
        value=True,
    )

    st.markdown("---")

    if st.button("▶ Predict Mutant Structure & Run Pipeline"):

        if not mutant_cam and not mutant_pep:
            st.error("❌ Please provide at least one mutant sequence.")
            st.stop()

        # Basic validation: only standard amino acids
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        for label, seq in [("CaM", mutant_cam), ("Peptide", mutant_pep)]:
            if seq:
                invalid = set(seq.upper()) - valid_aa
                if invalid:
                    st.error(
                        f"❌ Mutant {label} sequence contains invalid characters: "
                        f"{', '.join(sorted(invalid))}"
                    )
                    st.stop()

        # ── Step A: Predict mutant complex with Boltz-2 ───────────────────
        with st.spinner("Running Boltz-2 structure prediction (this may take 1–5 minutes)..."):
            try:
                mutant_pdb = logic5.predict_mutant_complex(
                    wildtype_pdb_path        = wt_pdb,
                    cam_chain_id             = cam_chain,
                    peptide_chain_id         = pep_chain,
                    mutant_cam_sequence      = mutant_cam     if mutant_cam     else None,
                    mutant_peptide_sequence  = mutant_pep     if mutant_pep     else None,
                    use_msa_server           = use_msa,
                )
                st.success(f"✅ Mutant structure predicted: `{mutant_pdb}`")
            except Exception as e:
                st.error(f"❌ Boltz-2 prediction failed: {e}")
                st.stop()

        # ── Step B: Fix protonation ───────────────────────────────────────
        with st.spinner("Fixing protonation..."):
            try:
                fixed_mutant_pdb = logic.fix_pdb_protonation(mutant_pdb)
                st.session_state["mutant_path"] = fixed_mutant_pdb
            except Exception as e:
                st.error(f"❌ Protonation fix failed: {e}")
                st.stop()

        # ── Step C: Run prppi on the mutant ──────────────────────────────
        with st.spinner("Running PRPPI on mutant structure..."):
            try:
                mutant_json = logic.run_prppi(fixed_mutant_pdb, cutoff=5.0)
                if mutant_json:
                    st.success(f"✅ PRPPI complete: `{mutant_json}`")
                    st.session_state["mutant_json_path"] = mutant_json
                    st.session_state["json_path"]        = mutant_json
                    st.session_state["path"]             = fixed_mutant_pdb

                    with st.expander("🔬 View Mutant JSON"):
                        with open(mutant_json) as f:
                            import json as jsonlib
                            st.json(jsonlib.load(f))

                    st.info(
                        "✅ Mutant JSON saved to session. "
                        "Navigate to **Step 2** to set residue ranges, "
                        "then **Step 3** to run OSPREY."
                    )
                else:
                    st.error("❌ PRPPI did not produce a JSON file.")
            except Exception as e:
                st.error(f"❌ PRPPI failed: {e}")


# ── Sidebar footer ────────────────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.write("**Software Team (Davis Lab)**")