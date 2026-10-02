"""
convert_to_amber.py — CHARMM -> AMBER topology conversion via ParmEd
=======================================================================
MMPBSA.py (AmberTools) expects Amber-format prmtop/inpcrd. ParmEd can
convert an already-parameterized CHARMM structure directly (energies
are force-field-agnostic once resolved to per-term numbers, so this
is a legitimate, commonly used conversion — NOT a re-parameterization).

Run this from inside the NAMD job directory (needs the .psf, the
solvated .pdb, and the same par_all36_*.prm / toppar_water_ions.str
files NAMD used).
"""
import parmed as pmd
from parmed.charmm import CharmmParameterSet

params = CharmmParameterSet(
    "par_all36_prot.prm", "par_all36_cgenff.prm", "toppar_water_ions.str"
)
struct = pmd.load_file("6Y4O_sidechain_fixed_fixed_pymol_redesigned_solvated.psf")
struct.load_parameters(params, copy=True)
struct.coordinates = pmd.load_file("6Y4O_sidechain_fixed_fixed_pymol_redesigned_solvated.pdb").coordinates

struct.save("complex.prmtop", format="amber", overwrite=True)
struct.save("complex.inpcrd", overwrite=True)
print("Wrote complex.prmtop / complex.inpcrd")
print("IMPORTANT: open complex.prmtop and confirm the water/ion residue")
print("names (expect WAT / Na+ / Cl- or similar) match the strip_mask")
print("used in mmpbsa.in — ParmEd's Amber writer may rename them.")
