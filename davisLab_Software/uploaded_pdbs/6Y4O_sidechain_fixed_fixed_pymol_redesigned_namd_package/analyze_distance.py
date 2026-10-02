"""
analyze_distance.py — Davis Lab CaM-RyR2 Portal
================================================
Computes N-domain / C-domain Ca distance from the PRODUCTION stage
trajectory only (mixing in heating/equilibration frames would bias
the result). Replicates Fig 2B of the manuscript.

Requirements: pip install MDAnalysis
"""
import os, sys
import MDAnalysis as mda
import numpy as np

PSF_FILE = "6Y4O_sidechain_fixed_fixed_pymol_redesigned_solvated.psf"
DCD_FILE = "04_production.dcd"
OUT_FILE = "domain_distance.dat"

for f in [PSF_FILE, DCD_FILE]:
    if not os.path.isfile(f):
        print(f"ERROR: Cannot find {f} — run the full staged protocol first.")
        sys.exit(1)

u = mda.Universe(PSF_FILE, DCD_FILE)
print(f"Frames: {len(u.trajectory)}  |  Atoms: {len(u.atoms)}")

sel_n = u.select_atoms("resid 40 and name CA and protein")
sel_c = u.select_atoms("resid 113 and name CA and protein")

frames, distances = [], []
for ts in u.trajectory:
    diff = sel_n.positions[0] - sel_c.positions[0]
    distances.append(float(np.linalg.norm(diff) / 10.0))
    frames.append(ts.frame)

with open(OUT_FILE, "w") as f:
    f.write("# Frame   Distance_nm\n")
    for fr, d in zip(frames, distances):
        f.write(f"{fr}   {d:.4f}\n")

arr = np.array(distances)
print(f"Saved {OUT_FILE}  |  mean={arr.mean():.3f} nm  min={arr.min():.3f}  max={arr.max():.3f}")
print("<1.0 nm ~ annealed | >1.5 nm ~ unlocked (RCaM1-like)")
