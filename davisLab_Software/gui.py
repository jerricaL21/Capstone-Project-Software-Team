import streamlit as st
import pandas as pd
import time
import os

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
    st.info("Upload the PDB file for the CaM-RyR2 complex.")
    
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
        st.subheader("PDB File Preview (Header)")
        st.code(pdb_bytes[:500] + "...", language="text")
        
        # Placeholder for 3D Visualization
        st.warning("Note: To view the 3D structure, you will need to install 'stmol'.")

        # ✅ User must click button — analysis only runs on explicit click
        if st.button("▶ Run Analysis"):
            import logic
            with st.spinner("Running PRPPI analysis..."):
                json_path = logic.run_prppi(path, cutoff=5.0)
            if json_path:
                st.success(f"Analysis complete! Results saved to: {json_path}")
                st.session_state["json_path"] = json_path  # save for Step 3

#Step 2: JSON Processing
elif workflow_step == "2. JSON Processing":
    st.header("Step 2: JSON Processing")

    #Checks if Step 1 already produced a JSON and saved it to memory. If not found, returns None.
    json_path = st.session_state.get("json_path", None)

    #If Step 1 ran successfully, show a green success message and wrap the single path in a list so the rest of the code can treat it the same as multiple files.
    if json_path:
        st.success(f"JSON ready from Step 1: {json_path}")
        json_paths = [json_path]
    else:
        uploaded_jsons = st.file_uploader( #Shows a file upload widget, restricts to JSON files only.Allows uploading more than one at once.
            "Upload Interaction JSON(s)",
            type=['json'],
            accept_multiple_files=True
        )
        json_paths = [] #Creates an empty list
        if uploaded_jsons: # Go ahead if user aploaded something
            os.makedirs("uploaded_jsons", exist_ok=True)# Creates folder called . Won't crash if folder exists
            for uploaded_json in uploaded_jsons:# Loops through each uploaded file one by one
                save_path = os.path.abspath(os.path.join("uploaded_jsons", uploaded_json.name))# File path of each json file
                with open(save_path, "wb") as f:# Saves the uploaded file to disk. "wb" means write in binary mode. .getvalue() gets the raw file bytes from the uploader.
                    f.write(uploaded_json.getvalue())
                json_paths.append(save_path)
            st.success(f"{len(json_paths)} JSON file(s) loaded:")
            for p in json_paths:
                st.write(f"• {p}")   
            st.session_state["json_paths"] = json_paths #Saves the list of paths to memory so other steps can access them.
    
    if json_paths:
        import logic2
        st.subheader("Set Residue Range Boundaries")
        mutant_start = st.text_input("Mutant Beginning")
        mutant_end   = st.text_input("Mutant Ending")
        ligand_start = st.text_input("Ligand Beginning")
        ligand_end   = st.text_input("Ligand Ending")

        if st.button("▶ Process JSON(s)"):
            for json_file in json_paths:#Loops through every JSON file and runs Part 1 on each one — adds Complex_Size metadata and sorts residues by distance.
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

    if st.button("Generate & Run OSPREY Script"):
        import logic3
        with st.spinner("Generating OSPREY Script..."):
            generated = logic3.create_folders_and_files(
                json_path, pdb_path, epsilon=epsilon, cpu_cores=int(cores)
            )
        st.success(f"{len(generated)} OSPREY script(s) generated!")
        for f in generated:
            st.write(f"• {f}")
        st.session_state["osprey_scripts"] = generated

#Step 4: Results
elif workflow_step == "4. Results & Analytics":
    st.header("Step 4: Redesign Results")
    
    # Mock data to show off the visualization
    data = {
        'Mutation': ['WT', 'ALA', 'VAL', 'ILE', 'LEU'],
        'K* Score': [10.2, 12.5, 9.8, 14.1, 11.2],
        'Stability': ['Stable', 'Stable', 'Unstable', 'Very Stable', 'Stable']
    }
    df = pd.DataFrame(data)
    df['Delta K*'] = df['K* Score'] - 10.2

    st.subheader("Binding Affinity Ranking (ΔK*)")
    st.bar_chart(df, x="Mutation", y="Delta K*")
    
    st.subheader("Sequence Details")
    st.table(df)

#Sidebar/Footer
st.sidebar.markdown("---")
st.sidebar.write("**Software Team (Davis Lab)**")