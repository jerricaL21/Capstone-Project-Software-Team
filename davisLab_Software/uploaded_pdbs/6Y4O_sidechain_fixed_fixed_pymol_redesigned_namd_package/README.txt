Davis Lab CaM-RyR2 Portal — staged NAMD/CHARMM36 package for 6Y4O_sidechain_fixed_fixed_pymol_redesigned
========================================================================

This reproduces (as closely as NAMD/CHARMM36 can reproduce a
GROMACS/AMBER99SB-ILDN study — see caveats in logic6.py's docstring)
the manuscript's explicit-solvent MD + standard equilibration
protocol.

BEFORE YOU RUN: download real CHARMM36 parameter files and place
them in this folder (these are the actual param/topology files —
the OpenMM Python force field used to BUILD the system is separate
from the files NAMD itself reads at runtime):
    par_all36_prot.prm
    par_all36_cgenff.prm
    toppar_water_ions.str
Source: http://mackerell.umaryland.edu/charmm_ff.shtml (registration
required) or via a CHARMM-GUI-generated toppar bundle.

Stages (run in order — each depends on the previous stage's restart
files, already wired up via bincoordinates/binvelocities/extendedSystem):
    01_minimize.namd    10,000-step restrained minimization
    02_heat.namd        NVT heating 0 -> 297.15 K, restrained
    03_equil.namd       NPT equilibration, restraint relaxed
    04_production.namd  NPT production, 100.0 ns, UNRESTRAINED

Run all four with:
    bash run_all_stages.sh
(edit the NAMD binary path/module load lines for your cluster first)

Then analyze:
    python analyze_distance.py

System summary
---------------
Atoms          : 32,155
Box (Angstrom) : 68.8 x 68.8 x 68.8
CMAP terms     : 158 (CHARMM36 backbone correction — should be > 0)
Protein (chain) residues, incl. any bound cofactor ions on the same chain letter: 1-142
Ligand/peptide residues: 143-166
NOTE: these are NOT necessarily one contiguous range — a bound cofactor
(e.g. Ca2+) can carry the protein's chain letter while sitting after the
ligand chain in residue-index order. The ranges above are exact.

Production length (100.0 ns) and temperature/pressure/ion
settings above match the manuscript's SI Appendix (MD Methods).
The SI does NOT state a replica count per system — run 2-3
independent replicates (different starting velocities) before
treating any single trajectory as representative.

For binding free energy (this package's MM/PBSA-equivalent step),
see mmpbsa/ subfolder and its README.
WARNINGS from system build: Detected incomplete exceptions. Not supported.; Unsupported Force type CustomNonbondedForce; Unsupported Force type CustomBondForce
