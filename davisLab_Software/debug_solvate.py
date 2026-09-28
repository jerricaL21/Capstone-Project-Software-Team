"""
debug_solvate.py — run this directly (NOT through the Streamlit GUI) to see
the real error that's happening inside solvate_and_build_system(), instead
of the generic "file not found" message the GUI shows.

Usage (from an Anaconda Prompt with pymol-env activated):

    python debug_solvate.py

Edit PDB_PATH and OUTPUT_DIR below to match your actual file locations
before running.
"""

import traceback
from logic7 import solvate_and_build_system

# --- EDIT THESE TWO PATHS to match your setup ------------------------------
PDB_PATH = r"D:\BME 4901.2 Capstone Project Software\Jerrica's Branch of Capstone-Project-Software-Team\Capstone-Project-Software-Team\davisLab_Software\uploaded_pdbs\6Y4O_sidechain_fixed_fixed_pymol_redesigned.pdb"
OUTPUT_DIR = r"D:\BME 4901.2 Capstone Project Software\Jerrica's Branch of Capstone-Project-Software-Team\Capstone-Project-Software-Team\davisLab_Software\uploaded_pdbs"
# ---------------------------------------------------------------------------

try:
    manifest = solvate_and_build_system(PDB_PATH, OUTPUT_DIR)
    print("SUCCESS!")
    print(manifest)
except Exception:
    print("=" * 70)
    print("REAL ERROR (this is what the GUI was hiding):")
    print("=" * 70)
    traceback.print_exc()