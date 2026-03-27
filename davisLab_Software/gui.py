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

importlib.reload(logic)
importlib.reload(logic2)
importlib.reload(logic3)
importlib.reload(logic4)

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
    ["1. Data Input", "2. JSON Processing", "3. OSPREY Execution", "4. Results & Analytics"])

#Step 1: Data Input
if workflow_step == "1. Data Input":
    st.header("Step 1: Structural Data Input")
    
    # Updated Uploader for PDB files
    uploaded_pdb = st.file_uploader("Upload Target PDB File", type=['pdb'])
    
    # Input for PDB ID (Optional backup)
    pdb_id = st.text_input("OR Enter Target PDB ID to fetch from RCSB", "2RL0")
    
    if uploaded_pdb is not None:
        
        # NEW: create folder if it doesn't exist
        os.makedirs("uploaded_pdbs", exist_ok=True)
        
        # NEW: define the path variable
        path = os.path.abspath(os.path.join("uploaded_pdbs", uploaded_pdb.name))
        
        # NEW: save the file to disk
        with open(path, "wb") as f:
            f.write(uploaded_pdb.getvalue())
        
        # Save path to session state so logic.py can access it
        st.session_state["path"] = path

        # Read the file content
        pdb_bytes = uploaded_pdb.getvalue().decode("utf-8")
        
        st.success(f"Successfully loaded: {uploaded_pdb.name}")
        
        # NEW: show saved path
        st.write(f"Saved file path: {path}")
        
        # Quick summary of the file
        with st.expander("📄 View PDB File Header"):
            st.code(pdb_bytes[:500] + "...", language="text")
        
        # Placeholder for 3D Visualization
        st.warning("Note: To view the 3D structure, you will need to install 'stmol'.")

        # ✅ User must click button — analysis only runs on explicit click
        if st.button("▶ Run Analysis"):
            with st.spinner("Fixing PDB protonation..."):
                fixed_pdb = logic.fix_pdb_protonation(path)
                st.session_state["path"] = fixed_pdb
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

#Step 2: JSON Processing
elif workflow_step == "2. JSON Processing":
    st.header("Step 2: JSON Processing")

    json_path = st.session_state.get("json_path", None)

    if json_path:
        st.success(f"JSON ready from Step 1: {json_path}")
        st.session_state["json_paths"] = [json_path]
        st.session_state["json_path"] = json_path  # re-save explicitly
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

        with st.form("residue_form"):
            mutant_start = st.text_input("Mutant Beginning")
            mutant_end   = st.text_input("Mutant Ending")
            ligand_start = st.text_input("Ligand Beginning")
            ligand_end   = st.text_input("Ligand Ending")
            submitted = st.form_submit_button("▶ Process JSON(s)")

        if submitted:
            if not all([mutant_start, mutant_end, ligand_start, ligand_end]):
                st.error("❌ Please fill in all four residue range fields before proceeding.")
            else:
                for json_file in json_paths:
                    logic2.add_variables_to_json(json_file, mutant_start, mutant_end, ligand_start, ligand_end)
                st.success("All JSON files processed! Proceed to Step 3.")


                
#Step 3: OSPREY Execution
elif workflow_step == "3. OSPREY Execution":
    st.header("Step 3: OSPREY Simulation Setup")
    
    json_path = st.session_state.get("json_path", None)
    pdb_path  = st.session_state.get("path", None)

    if not json_path or not pdb_path:
        st.warning("Please complete Steps 1 and 2 first.")
        st.stop()

    #Load Residues from JSON
    import json as jsonlib
    with open(json_path) as f:
        data = jsonlib.load(f)
    residues = [k for k in data.keys() if k != "Complex_Size"]
        
    
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Design Site")
        target_res = st.selectbox("Select Residue to Mutate", residues)
        mutations = st.multiselect("Select Amino Acids", 
        ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'GLU', 'TYR', 'ASP'])
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
#Sidebar/Footer
st.sidebar.markdown("---")
st.sidebar.write("**Software Team (Davis Lab)**")