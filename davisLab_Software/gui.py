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
        # Get the directory where results are stored
        base_dir = os.path.dirname(os.path.dirname(osprey_scripts[0]))

        # Run the updated analysis logic
        pivot_delta, summary_table, error = logic4.analyze_results(base_dir)

        if error:
            st.error(f"Analysis Error: {error}")
        else:
            # --- 1. Top Level Metrics ---
            st.subheader("📊 Affinity Summary")
            
            # Use columns to show the "Best" result at a glance
            col1, col2 = st.columns([2, 1])
            
            with col1:
                st.markdown("**Comparison of Mutations vs Wild-Type ($\Delta K^*$)**")
                # Grouped bar chart
                # Since pivot_delta is rounded in logic4, these bars look clean
                st.bar_chart(pivot_delta)
            
            with col2:
                # Find the sequence with the highest K* improvement (Max Delta)
                try:
                    best_val = pivot_delta.max().max()
                    best_mut = pivot_delta.max().idxmax()
                    st.metric("Top Improvement", f"+{best_val:.2f}", f"Amino Acid: {best_mut}")
                except:
                    st.info("Insufficient data for metrics")

            st.divider()

            # --- 2. Detailed Data Tables ---
            st.subheader("🔬 Detailed Sequence Scores")
            
            tab1, tab2 = st.tabs(["Summary Table", "Raw Delta Matrix"])
            
            with tab1:
                st.write("Overview of all tested mutations and their absolute $K^*$ scores:")
                # Display the clean summary table we built
                st.dataframe(summary_table, use_container_width=True, hide_index=True)
            
            with tab2:
                st.write("Matrix of $\Delta K^*$ scores (0.00 = Wild Type reference):")
                # Display the pivot table with highlighting
                st.dataframe(pivot_delta.style.background_gradient(cmap='RdYlGn'), use_container_width=True)

            # --- 3. Export Section ---
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