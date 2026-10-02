Binding free energy — MM/PBSA-equivalent step
===============================================
This is the NAMD/AMBER-world analogue of the manuscript's
GROMACS g_mmpbsa calculation. It requires AmberTools (MMPBSA.py)
installed on whatever machine you run it on — not included here.

Masks (computed from your structure's chain IDs):
    receptor (CaM + bound cofactor ions): :1-142
    ligand   (RyR2 peptide):                                   :143-166

NOTE: MMPBSA.py's single-trajectory method also needs separate
receptor.prmtop / ligand.prmtop files (same atoms, subset topology).
Generate those with ParmEd's `struct[receptor_mask]` / `struct[ligand_mask]`
slicing in convert_to_amber.py, or via `ante-MMPBSA.py`, before running
run_mmpbsa.sh — the driver script assumes they already exist.

Steps:
    1. cd into the mmpbsa/ folder (this folder)
    2. Copy or symlink the .psf, solvated .pdb, and CHARMM36 .prm/.str
       files from the NAMD job directory into here
    3. bash run_mmpbsa.sh
