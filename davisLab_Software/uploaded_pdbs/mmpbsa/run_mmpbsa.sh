#!/bin/bash
# Requires AmberTools (MMPBSA.py) installed on this machine.
set -e
python convert_to_amber.py
MMPBSA.py -O -i mmpbsa.in \
    -cp complex.prmtop \
    -rp receptor.prmtop -lp ligand.prmtop \
    -y ../6Y4O_sidechain_fixed_fixed_pymol_redesigned_namd_package/04_production.dcd
echo "See FINAL_RESULTS_MMPBSA.dat for dG_bind"
