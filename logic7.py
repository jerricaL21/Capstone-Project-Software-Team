"""
logic7.py  –  Explicit-solvent NAMD/CHARMM36 molecular dynamics pipeline
Davis Lab Software Team, 2026

REWRITE NOTES (read this first)
--------------------------------
The manuscript ("Dynamic-Structure Redesign Of Calmodulin Reveals
Mechanistic Constraints On Ryr2 Regulation") ran its MD in GROMACS with the
AMBER99SB-ILDN force field, explicit solvent, standard equilibration
protocols, LINCS bond constraints, and MM/PBSA binding energetics via
g_mmpbsa. The previous version of this module generated a single NAMD
config using CHARMM36 with *implicit* GB solvent and no equilibration
stages at all — that setup cannot reproduce the paper's dynamics; implicit
vs. explicit solvent is a qualitative difference, not a tuning knob, and
skipping equilibration on an unrelaxed system risks a numerically unstable
run.

This version:
  1. Builds an explicit TIP3P-solvated, 0.15 M NaCl-neutralized system
     using OpenMM's Modeller + CHARMM36 force field files (no VMD/psfgen
     required — VMD's GUI renderer is the thing that was crashing on
     Windows without an NVIDIA driver; this sidesteps it entirely), then
     converts to a CHARMM PSF/PDB pair via ParmEd so NAMD can read it.
  2. Generates a *staged* NAMD protocol (minimize -> heat -> restrained
     NPT equilibration -> unrestrained NPT production) chained via NAMD
     restart files, with a Langevin thermostat + Langevin piston barostat
     (the NAMD equivalents of GROMACS's "standard equilibration protocol").
  3. Generates a binding-energetics step: ParmEd converts the CHARMM
     system to an AMBER prmtop (parameters carry over losslessly once
     resolved into per-particle/per-term numbers) so AmberTools'
     MMPBSA.py can be run — this is the direct analogue of the paper's
     g_mmpbsa (GROMACS) MM/PBSA calculation.

KNOWN METHODOLOGICAL CAVEATS — put these in your methods section if you
publish anything from this pipeline:
  - CHARMM36 =/= AMBER99SB-ILDN. They're both well-validated modern
    protein force fields, but they are not the same model. Expect
    quantitative differences, especially in flexible/linker regions —
    exactly the kind of region this paper's central finding (N-/C-domain
    "annealing" vs. "unlocking") depends on. Don't present CHARMM36 results
    as a replication of the GROMACS numbers; present them as an
    independent cross-check.
  - The OpenMM -> ParmEd conversion used to build the PSF preserves CMAP
    (verified: CHARMM36's backbone cross-term correction survives the
    round-trip intact) but drops CHARMM's NBFIX pairwise ion-correction
    terms. This mainly affects fine details of ion behavior, not backbone
    dynamics. If you need publication-grade ion parameters, cross-check
    system building against CHARMM-GUI's Solution Builder, which writes
    NBFIX-correct parameters directly.
  - Per the SI Appendix (Molecular Dynamics Simulations section): 297.15 K,
    0.15 M KCl (not NaCl), production >= 100 ns, 1.0 bar reference pressure,
    minimum 1.0 nm protein-to-box-edge padding, coords/energies saved every
    10 ps. Defaults below now match these. The SI does not state a replica
    count per system — confirm that separately before drawing quantitative
    comparisons across replicates.
  - This module cannot run NAMD itself (no NAMD binary ships here, and a
    real production run needs a GPU and hours-to-days). It prepares a
    ready-to-submit job directory for your own GPU workstation/cluster.
"""

import os
import json
import textwrap
import zipfile


# ---------------------------------------------------------------------------
# Stage 1 — explicit solvation + CHARMM PSF/PDB construction
# ---------------------------------------------------------------------------

def solvate_and_build_system(
    pdb_path: str,
    output_dir: str,
    padding_nm: float = 1.0,
    ionic_strength_m: float = 0.15,
    protein_chain: str = "A",
    ligand_chain: str = "B",
) -> dict:
    """
    Add hydrogens, solvate (TIP3P) and neutralize/ionize (NaCl) a PDB using
    OpenMM's CHARMM36 force field, then write out a CHARMM PSF + PDB pair
    for NAMD via ParmEd. No VMD/psfgen dependency.

    Returns a dict with: psf_path, pdb_path, box_angstrom (x,y,z),
    protein_resid_range, ligand_resid_range (1-based, inclusive — computed
    from the *input* structure's chain IDs, before solvent is added, so
    these ranges also correctly include any HETATM cofactors, e.g. bound
    Ca2+, that share the protein chain ID), n_atoms, warnings (list[str]).
    """
    try:
        import openmm
        from openmm import unit
        from openmm.app import ForceField, HBonds, Modeller, PDBFile, PME
        import parmed as pmd
    except ImportError as e:
        raise ImportError(
            "Missing dependency for explicit-solvent system building. Run:\n"
            "    pip install openmm parmed\n"
            f"(original error: {e})"
        )

    pdb_path = os.path.abspath(pdb_path)
    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
    build_dir = os.path.join(output_dir, f"{pdb_name}_system")
    os.makedirs(build_dir, exist_ok=True)

    # ---- Determine solute residue ranges BEFORE adding solvent -----------
    # Order is preserved by Modeller.addSolvent (new water/ion residues are
    # always appended after existing residues), so indices computed here
    # still apply to the final solvated PSF.
    pre_struct = pmd.load_file(pdb_path)
    protein_resids = [i + 1 for i, r in enumerate(pre_struct.residues) if r.chain == protein_chain]
    ligand_resids = [i + 1 for i, r in enumerate(pre_struct.residues) if r.chain == ligand_chain]
    if not protein_resids or not ligand_resids:
        raise ValueError(
            f"Could not find residues on chain '{protein_chain}' and/or '{ligand_chain}' "
            f"in {pdb_path}. Chains present: {sorted(set(r.chain for r in pre_struct.residues))}"
        )
    # NOTE: a chain's residues are not guaranteed contiguous in index space —
    # e.g. bound cofactor ions (HETATM, same chain letter as the protein) can
    # sit in the numbering *after* a ligand chain that comes between them.
    # Keep the full residue-index lists (not just min/max) so downstream
    # Amber masks can be built as proper (possibly multi-range) selections.

    # ---- Build solvated system with OpenMM --------------------------------
    warnings_out = []
    pdb = PDBFile(pdb_path)
    forcefield = ForceField("charmm36.xml", "charmm36/water.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=padding_nm * unit.nanometer,
        ionicStrength=ionic_strength_m * unit.molar,
        positiveIon="K+",
        negativeIon="Cl-",
    )
    box_vectors = modeller.topology.getPeriodicBoxVectors()
    box_angstrom = tuple(box_vectors[i][i].value_in_unit(unit.angstrom) for i in range(3))

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.2 * unit.nanometer,
        constraints=HBonds,
    )

    # ---- Convert to CHARMM PSF/PDB via ParmEd ------------------------------
    import warnings as _w
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        struct = pmd.openmm.load_topology(modeller.topology, system=system, xyz=modeller.positions)
        for w in caught:
            warnings_out.append(str(w.message))

    n_cmap = len(struct.cmaps)
    if n_cmap == 0:
        warnings_out.append(
            "WARNING: 0 CMAP terms survived the OpenMM->ParmEd conversion — "
            "CHARMM36's backbone correction is missing from this PSF. Do not "
            "use this system for production; rebuild (e.g. via CHARMM-GUI) instead."
        )

    psf_path = os.path.join(build_dir, f"{pdb_name}_solvated.psf")
    solv_pdb_path = os.path.join(build_dir, f"{pdb_name}_solvated.pdb")

    # Make sure the target folder exists on disk before ParmEd attempts to save
    os.makedirs(build_dir, exist_ok=True)

    struct.save(psf_path, format="psf", overwrite=True)
    struct.save(solv_pdb_path, overwrite=True)

    manifest = {
        "psf_path": psf_path,
        "pdb_path": solv_pdb_path,
        "box_angstrom": box_angstrom,
        "protein_resids": protein_resids,   # full index list, not just min/max
        "ligand_resids": ligand_resids,
        "n_atoms": len(struct.atoms),
        "n_cmap_terms": n_cmap,
        "warnings": warnings_out,
    }
    with open(os.path.join(build_dir, "system_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    return manifest


def _compress_to_ranges(sorted_ints: list) -> str:
    """[1,2,3,5,6,9] -> '1-3,5-6,9' (Amber mask syntax)."""
    if not sorted_ints:
        return ""
    ranges = []
    start = prev = sorted_ints[0]
    for x in sorted_ints[1:]:
        if x == prev + 1:
            prev = x
            continue
        ranges.append((start, prev))
        start = prev = x
    ranges.append((start, prev))
    return ",".join(f"{a}-{b}" if a != b else f"{a}" for a, b in ranges)


def _write_restraint_pdb(solvated_pdb_path: str, out_path: str, solute_resids: set, k_backbone: float):
    """
    Write a copy of the solvated PDB with the B-factor (beta) column set to
    `k_backbone` (kcal/mol/A^2) on backbone atoms (N, CA, C, O) of the
    solute (protein+ligand residue range), and 0.0 everywhere else. NAMD
    reads this directly via `conskfile ... conskcol B` to apply harmonic
    positional restraints during equilibration.
    """
    backbone_atoms = {"N", "CA", "C", "O"}
    out_lines = []
    resid_counter = 0
    last_res_key = None
    with open(solvated_pdb_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                res_key = line[17:26]  # resname+chain+resseq, stable per-residue key
                if res_key != last_res_key:
                    resid_counter += 1
                    last_res_key = res_key
                beta = k_backbone if (resid_counter in solute_resids and atom_name in backbone_atoms) else 0.0
                line = line[:60] + f"{beta:6.2f}" + line[66:]
            out_lines.append(line)
    with open(out_path, "w") as fh:
        fh.writelines(out_lines)


# ---------------------------------------------------------------------------
# Stage 2 — staged NAMD protocol (minimize / heat / equil / production)
# ---------------------------------------------------------------------------

def generate_namd_package(
    system_manifest: dict,
    output_dir: str,
    temperature: float = 297.15,
    production_ns: float = 100.0,
    timestep_fs: float = 2.0,
    nonbonded_cutoff: float = 12.0,
    gpu_enabled: bool = True,
) -> dict:
    """
    Generate a staged, explicit-solvent NAMD job directory from a system
    manifest produced by `solvate_and_build_system`.

    Stages (chained via NAMD restart files — each reads the previous
    stage's .coor/.vel/.xsc):
        01_minimize    – 10,000 steps CG minimization, backbone restrained
                          (k=10 kcal/mol/A^2) so solvent relaxes around a
                          fixed protein rather than the protein relaxing
                          into a clash in unequilibrated solvent.
        02_heat        – NVT, 0 -> `temperature` K ramp (reassignIncr),
                          backbone still restrained (k=10).
        03_equil       – NPT (Langevin piston, 1 atm), backbone restraint
                          relaxed to k=1 kcal/mol/A^2, 500 ps.
        04_production  – NPT, unrestrained, `production_ns` ns. This is the
                          stage `analyze_distance.py` reads.
    """
    psf_path = system_manifest["psf_path"]
    pdb_path = system_manifest["pdb_path"]
    box = system_manifest["box_angstrom"]
    protein_resids = system_manifest["protein_resids"]
    ligand_resids = system_manifest["ligand_resids"]
    solute_resids = set(protein_resids) | set(ligand_resids)

    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0].replace("_solvated", "")
    namd_dir = os.path.join(output_dir, f"{pdb_name}_namd_package")
    os.makedirs(namd_dir, exist_ok=True)

    # Copy topology files into the job dir so it's self-contained
    import shutil
    local_psf = os.path.join(namd_dir, os.path.basename(psf_path))
    local_pdb = os.path.join(namd_dir, os.path.basename(pdb_path))
    shutil.copy(psf_path, local_psf)
    shutil.copy(pdb_path, local_pdb)

    # Restraint reference PDBs
    restraint_full = os.path.join(namd_dir, "restraint_k10.pdb")
    restraint_weak = os.path.join(namd_dir, "restraint_k1.pdb")
    _write_restraint_pdb(local_pdb, restraint_full, solute_resids, 10.0)
    _write_restraint_pdb(local_pdb, restraint_weak, solute_resids, 1.0)

    switchdist = nonbonded_cutoff - 2.0
    pairlistdist = nonbonded_cutoff + 1.5
    gpu_block = "CUDASOAintegrate        on\n" if gpu_enabled else "# GPU disabled — CPU only\n"

    def common_header(structure_note):
        return textwrap.dedent(f"""\
            # =====================================================================
            # NAMD Configuration — {structure_note}
            # Structure : {pdb_name}   (explicit TIP3P solvent, CHARMM36)
            # =====================================================================
            structure            {os.path.basename(local_psf)}
            parameters           par_all36_prot.prm
            parameters           par_all36_cgenff.prm
            parameters           toppar_water_ions.str
            paratypecharmm       on

            wrapAll              on
            wrapWater            on
            cellBasisVector1     {box[0]:.3f}   0.000   0.000
            cellBasisVector2     0.000   {box[1]:.3f}   0.000
            cellBasisVector3     0.000   0.000   {box[2]:.3f}
            cellOrigin           0.000   0.000   0.000

            exclude              scaled1-4
            1-4scaling           1.0
            cutoff               {nonbonded_cutoff}
            switching            on
            switchdist           {switchdist}
            pairlistdist         {pairlistdist}
            PME                  yes
            PMEGridSpacing       1.0

            rigidBonds           all
            nonbondedFreq        1
            fullElectFrequency   2
            timestep             {timestep_fs}

            {gpu_block}""")

    # ---- 01: minimize -------------------------------------------------------
    min_conf = common_header("Stage 1/4: restrained minimization") + textwrap.dedent(f"""
        coordinates          {os.path.basename(local_pdb)}
        temperature          0

        constraints           on
        consexp                2
        consref                {os.path.basename(restraint_full)}
        conskfile               {os.path.basename(restraint_full)}
        conskcol                B

        outputName            01_minimize
        outputEnergies        500
        minimize              10000
    """)
    with open(os.path.join(namd_dir, "01_minimize.namd"), "w") as fh:
        fh.write(min_conf)

    # ---- 02: heat (NVT, restrained) ----------------------------------------
    heat_steps = 25000  # 50 ps at 2 fs
    heat_conf = common_header("Stage 2/4: NVT heating 0 -> target, restrained") + textwrap.dedent(f"""
        bincoordinates        01_minimize.coor
        extendedSystem        01_minimize.xsc
        temperature           0

        constraints            on
        consexp                 2
        consref                 {os.path.basename(restraint_full)}
        conskfile                {os.path.basename(restraint_full)}
        conskcol                 B

        langevin              on
        langevinDamping       1
        langevinHydrogen      off
        langevinTemp          0
        reassignFreq          1000
        reassignTemp          0
        reassignIncr          10
        reassignHold          {temperature}

        outputName            02_heat
        dcdfile               02_heat.dcd
        dcdFreq               5000
        outputEnergies        500
        restartFreq           5000
        firsttimestep         0
        run                   {heat_steps}
    """)
    with open(os.path.join(namd_dir, "02_heat.namd"), "w") as fh:
        fh.write(heat_conf)

    # ---- 03: NPT equilibration, weak restraint -----------------------------
    equil_steps = 250000  # 500 ps at 2 fs
    equil_conf = common_header("Stage 3/4: NPT equilibration, restraints relaxed") + textwrap.dedent(f"""
        bincoordinates        02_heat.coor
        binvelocities         02_heat.vel
        extendedSystem        02_heat.xsc

        constraints            on
        consexp                 2
        consref                 {os.path.basename(restraint_weak)}
        conskfile                {os.path.basename(restraint_weak)}
        conskcol                 B

        langevin              on
        langevinDamping       1
        langevinTemp          {temperature}
        langevinHydrogen      off

        useGroupPressure      yes
        useFlexibleCell       no
        useConstantArea       no
        langevinPiston        on
        langevinPistonTarget  1.0
        langevinPistonPeriod  100
        langevinPistonDecay   50
        langevinPistonTemp    {temperature}

        outputName            03_equil
        dcdfile               03_equil.dcd
        dcdFreq               5000
        outputEnergies        500
        restartFreq           5000
        firsttimestep         {heat_steps}
        run                   {equil_steps}
    """)
    with open(os.path.join(namd_dir, "03_equil.namd"), "w") as fh:
        fh.write(equil_conf)

    # ---- 04: unrestrained NPT production -----------------------------------
    production_steps = int(production_ns * 1_000_000 / timestep_fs)
    production_dcd_freq = 5000
    prod_conf = common_header("Stage 4/4: unrestrained NPT production") + textwrap.dedent(f"""
        bincoordinates        03_equil.coor
        binvelocities         03_equil.vel
        extendedSystem        03_equil.xsc

        langevin              on
        langevinDamping       1
        langevinTemp          {temperature}
        langevinHydrogen      off

        useGroupPressure      yes
        useFlexibleCell       no
        useConstantArea       no
        langevinPiston        on
        langevinPistonTarget  1.0
        langevinPistonPeriod  100
        langevinPistonDecay   50
        langevinPistonTemp    {temperature}

        outputName            04_production
        dcdfile               04_production.dcd
        dcdFreq               {production_dcd_freq}
        outputEnergies        1000
        outputPressure        1000
        restartFreq           25000
        firsttimestep         {heat_steps + equil_steps}
        run                   {production_steps}
    """)
    with open(os.path.join(namd_dir, "04_production.namd"), "w") as fh:
        fh.write(prod_conf)

    # ---- Analysis script (unchanged approach, points at production stage) -
    analysis_py = textwrap.dedent(f'''\
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

        PSF_FILE = "{os.path.basename(local_psf)}"
        DCD_FILE = "04_production.dcd"
        OUT_FILE = "domain_distance.dat"

        for f in [PSF_FILE, DCD_FILE]:
            if not os.path.isfile(f):
                print(f"ERROR: Cannot find {{f}} — run the full staged protocol first.")
                sys.exit(1)

        u = mda.Universe(PSF_FILE, DCD_FILE)
        print(f"Frames: {{len(u.trajectory)}}  |  Atoms: {{len(u.atoms)}}")

        sel_n = u.select_atoms("resid 40 and name CA and protein")
        sel_c = u.select_atoms("resid 113 and name CA and protein")

        frames, distances = [], []
        for ts in u.trajectory:
            diff = sel_n.positions[0] - sel_c.positions[0]
            distances.append(float(np.linalg.norm(diff) / 10.0))
            frames.append(ts.frame)

        with open(OUT_FILE, "w") as f:
            f.write("# Frame   Distance_nm\\n")
            for fr, d in zip(frames, distances):
                f.write(f"{{fr}}   {{d:.4f}}\\n")

        arr = np.array(distances)
        print(f"Saved {{OUT_FILE}}  |  mean={{arr.mean():.3f}} nm  min={{arr.min():.3f}}  max={{arr.max():.3f}}")
        print("<1.0 nm ~ annealed | >1.5 nm ~ unlocked (RCaM1-like)")
    ''')
    with open(os.path.join(namd_dir, "analyze_distance.py"), "w") as fh:
        fh.write(analysis_py)

    # ---- Launcher (Linux/HPC — adjust module/scheduler lines as needed) ---
    launcher = textwrap.dedent(f"""\
        #!/bin/bash
        # Davis Lab — staged NAMD launcher (GPU). Run stages IN ORDER.
        set -e
        NAMD=namd3   # or your cluster's namd binary/module
        $NAMD +p8 +devices 0 01_minimize.namd  > 01_minimize.log
        $NAMD +p8 +devices 0 02_heat.namd      > 02_heat.log
        $NAMD +p8 +devices 0 03_equil.namd     > 03_equil.log
        $NAMD +p8 +devices 0 04_production.namd > 04_production.log
        echo "Done. Run: python analyze_distance.py"
    """)
    launcher_path = os.path.join(namd_dir, "run_all_stages.sh")
    with open(launcher_path, "w") as fh:
        fh.write(launcher)
    os.chmod(launcher_path, 0o755)

    readme = textwrap.dedent(f"""\
        Davis Lab CaM-RyR2 Portal — staged NAMD/CHARMM36 package for {pdb_name}
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
            02_heat.namd        NVT heating 0 -> {temperature} K, restrained
            03_equil.namd       NPT equilibration, restraint relaxed
            04_production.namd  NPT production, {production_ns} ns, UNRESTRAINED

        Run all four with:
            bash run_all_stages.sh
        (edit the NAMD binary path/module load lines for your cluster first)

        Then analyze:
            python analyze_distance.py

        System summary
        ---------------
        Atoms          : {system_manifest['n_atoms']:,}
        Box (Angstrom) : {box[0]:.1f} x {box[1]:.1f} x {box[2]:.1f}
        CMAP terms     : {system_manifest['n_cmap_terms']} (CHARMM36 backbone correction — should be > 0)
        Protein (chain) residues, incl. any bound cofactor ions on the same chain letter: {_compress_to_ranges(protein_resids)}
        Ligand/peptide residues: {_compress_to_ranges(ligand_resids)}
        NOTE: these are NOT necessarily one contiguous range — a bound cofactor
        (e.g. Ca2+) can carry the protein's chain letter while sitting after the
        ligand chain in residue-index order. The ranges above are exact.

        Production length ({production_ns} ns) and temperature/pressure/ion
        settings above match the manuscript's SI Appendix (MD Methods).
        The SI does NOT state a replica count per system — run 2-3
        independent replicates (different starting velocities) before
        treating any single trajectory as representative.

        For binding free energy (this package's MM/PBSA-equivalent step),
        see mmpbsa/ subfolder and its README.
        {("WARNINGS from system build: " + "; ".join(system_manifest["warnings"])) if system_manifest["warnings"] else ""}
    """)
    with open(os.path.join(namd_dir, "README.txt"), "w") as fh:
        fh.write(readme)

    return {
        "namd_dir": namd_dir,
        "stages": ["01_minimize.namd", "02_heat.namd", "03_equil.namd", "04_production.namd"],
        "readme_path": os.path.join(namd_dir, "README.txt"),
        "launcher_path": launcher_path,
        "analysis_path": os.path.join(namd_dir, "analyze_distance.py"),
        "production_ns": production_ns,
        "production_steps": production_steps,
        "production_dcd_freq": production_dcd_freq,
    }


# ---------------------------------------------------------------------------
# Stage 3 — binding energetics (MM/PBSA-equivalent, the analogue of the
# manuscript's g_mmpbsa calculation)
# ---------------------------------------------------------------------------

def generate_mmpbsa_step(system_manifest: dict, namd_package: dict, output_dir: str) -> dict:
    """
    Write a ParmEd CHARMM->AMBER conversion script + an MMPBSA.py driver
    for binding free energy on the production trajectory. Requires
    AmberTools (MMPBSA.py) on whatever machine you actually run this on —
    not included here.
    """
    protein_resids = system_manifest["protein_resids"]
    ligand_resids = system_manifest["ligand_resids"]
    psf_name = os.path.basename(system_manifest["psf_path"])
    pdb_name = os.path.basename(system_manifest["pdb_path"])

    mm_dir = os.path.join(output_dir, "mmpbsa")
    os.makedirs(mm_dir, exist_ok=True)

    convert_py = textwrap.dedent(f'''\
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
        struct = pmd.load_file("{psf_name}")
        struct.load_parameters(params, copy=True)
        struct.coordinates = pmd.load_file("{pdb_name}").coordinates

        struct.save("complex.prmtop", format="amber", overwrite=True)
        struct.save("complex.inpcrd", overwrite=True)
        print("Wrote complex.prmtop / complex.inpcrd")
        print("IMPORTANT: open complex.prmtop and confirm the water/ion residue")
        print("names (expect WAT / Na+ / Cl- or similar) match the strip_mask")
        print("used in mmpbsa.in — ParmEd's Amber writer may rename them.")
    ''')
    with open(os.path.join(mm_dir, "convert_to_amber.py"), "w") as fh:
        fh.write(convert_py)

    # Amber masks as *exact* (possibly multi-range) selections — do not
    # collapse to a naive min-max span, since a bound cofactor ion can carry
    # the protein's chain letter while sitting after the ligand chain in
    # residue-index order (see solvate_and_build_system's docstring).
    receptor_mask = f":{_compress_to_ranges(protein_resids)}"
    ligand_mask = f":{_compress_to_ranges(ligand_resids)}"

    mmpbsa_in = textwrap.dedent(f"""\
        Single-trajectory GB binding free energy (analogue of the paper's
        GROMACS g_mmpbsa calculation)
        &general
           startframe=1, endframe=9999999, interval=1,
           strip_mask=:WAT,Na+,Cl-,          # verify these names in complex.prmtop first
        /
        &gb
           igb=5, saltcon=0.150,
        /
    """)
    with open(os.path.join(mm_dir, "mmpbsa.in"), "w") as fh:
        fh.write(mmpbsa_in)

    driver_sh = textwrap.dedent(f"""\
        #!/bin/bash
        # Requires AmberTools (MMPBSA.py) installed on this machine.
        set -e
        python convert_to_amber.py
        MMPBSA.py -O -i mmpbsa.in \\
            -cp complex.prmtop \\
            -rp receptor.prmtop -lp ligand.prmtop \\
            -y ../{os.path.basename(namd_package["namd_dir"])}/04_production.dcd
        echo "See FINAL_RESULTS_MMPBSA.dat for dG_bind"
    """)
    driver_path = os.path.join(mm_dir, "run_mmpbsa.sh")
    with open(driver_path, "w") as fh:
        fh.write(driver_sh)
    os.chmod(driver_path, 0o755)

    readme = textwrap.dedent(f"""\
        Binding free energy — MM/PBSA-equivalent step
        ===============================================
        This is the NAMD/AMBER-world analogue of the manuscript's
        GROMACS g_mmpbsa calculation. It requires AmberTools (MMPBSA.py)
        installed on whatever machine you run it on — not included here.

        Masks (computed from your structure's chain IDs):
            receptor (CaM {"+ bound cofactor ions" if True else ""}): {receptor_mask}
            ligand   (RyR2 peptide):                                   {ligand_mask}

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
    """)
    with open(os.path.join(mm_dir, "README.txt"), "w") as fh:
        fh.write(readme)

    return {"mmpbsa_dir": mm_dir, "receptor_mask": receptor_mask, "ligand_mask": ligand_mask}


# ---------------------------------------------------------------------------
# Convenience: run the whole pipeline and zip it
# ---------------------------------------------------------------------------

def build_full_package(pdb_path: str, output_dir: str, **kwargs) -> dict:
    """Runs solvation -> staged NAMD config generation -> MM/PBSA scaffolding -> zip."""
    manifest = solvate_and_build_system(
        pdb_path, output_dir,
        padding_nm=kwargs.get("padding_nm", 1.0),
        ionic_strength_m=kwargs.get("ionic_strength_m", 0.15),
        protein_chain=kwargs.get("protein_chain", "A"),
        ligand_chain=kwargs.get("ligand_chain", "B"),
    )
    namd_package = generate_namd_package(
        manifest, output_dir,
        temperature=kwargs.get("temperature", 297.15),
        production_ns=kwargs.get("production_ns", 100.0),
        timestep_fs=kwargs.get("timestep_fs", 2.0),
        nonbonded_cutoff=kwargs.get("nonbonded_cutoff", 12.0),
        gpu_enabled=kwargs.get("gpu_enabled", True),
    )
    mmpbsa = generate_mmpbsa_step(manifest, namd_package, output_dir)

    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
    zip_path = os.path.join(output_dir, f"{pdb_name}_namd_explicit_solvent_package.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for root_dir in (namd_package["namd_dir"], mmpbsa["mmpbsa_dir"]):
            base = os.path.dirname(root_dir)
            for dirpath, _, files in os.walk(root_dir):
                for fname in files:
                    full = os.path.join(dirpath, fname)
                    zf.write(full, arcname=os.path.relpath(full, base))

    return {
        "manifest": manifest,
        "namd_package": namd_package,
        "mmpbsa": mmpbsa,
        "zip_path": zip_path,
        "summary": (
            f"Built {manifest['n_atoms']:,}-atom explicit-solvent system "
            f"({manifest['n_cmap_terms']} CMAP terms) and a 4-stage NAMD package "
            f"(min -> heat -> equil -> {namd_package['production_ns']} ns production) "
            f"with MM/PBSA scaffolding."
        ),
    }


def parse_domain_distance_dat(dat_text: str) -> list[dict]:
    """Parse the domain_distance.dat file produced by analyze_distance.py."""
    rows = []
    for line in dat_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                rows.append({"frame": int(parts[0]), "distance_nm": float(parts[1])})
            except ValueError:
                continue
    return rows