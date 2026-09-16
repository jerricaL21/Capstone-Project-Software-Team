

import osprey
import os
import shutil
import glob

def _find_system_jvm_dll():
    """
    Try to locate jvm.dll from a system-installed JDK/JRE.
    Returns None if nothing suitable is found, in which case OSPREY
    falls back to its own bundled JRE.

    This is used because OSPREY's own bundled JRE can crash on some
    newer CPUs (a known old-JVM-build compatibility issue). If the
    user has ANY separate JDK installed (via JAVA_HOME, PATH, or a
    common install location), we prefer that instead.
    """
    candidates = []

    java_home = os.environ.get("JAVA_HOME")
    if java_home:
        candidates.append(java_home)

    java_exe = shutil.which("java")
    if java_exe:
        bin_dir = os.path.dirname(java_exe)
        candidates.append(os.path.dirname(bin_dir))

    candidates.extend(glob.glob(r"C:\Program Files\Java\jdk-*"))
    candidates.extend(glob.glob(r"C:\Program Files\Eclipse Adoptium\jdk-*"))

    for root in candidates:
        for rel in ("bin/server/jvm.dll", "bin/client/jvm.dll"):
            path = os.path.join(root, *rel.split("/"))
            if os.path.isfile(path):
                return path

    return None


# Optional manual override: set an OSPREY_JRE_PATH environment variable
# to point at a specific jvm.dll, otherwise auto-detect a system JDK.
_jvm_path = os.environ.get("OSPREY_JRE_PATH") or _find_system_jvm_dll()

if _jvm_path:
    # osprey.start() would normally pick OSPREY's own bundled JRE, which can
    # crash on some newer CPUs. Instead, call OSPREY's internal common-setup
    # routine directly, overriding which JRE it boots, while still getting all
    # the same setup osprey.start() does (classpath, the Java class factory
    # 'c', WILD_TYPE, Forcefield, etc.) - just pointed at a working JDK.
    osprey._start_jvm_common(lambda _ignored_default_path: osprey.jvm.start(
        _jvm_path,
        heapSizeMiB=4096,
        enableAssertions=False,
        stackSizeMiB=16,
        garbageSizeMiB=512,
        allowRemoteManagement=False,
        attachJvmDebugger=False,
    ))
    if osprey._print_preamble:
        print("Using up to 4096 MiB heap memory: 512 MiB for garbage, %d MiB for storage" % (4096 - 512))
else:
    # No system JDK found on this machine - fall back to OSPREY's bundled JRE
    osprey.start(heapSizeMiB=4096, garbageSizeMiB=512)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb("D:/BME 4901.2 Capstone Project Software/Jerrica's Branch of Capstone-Project-Software-Team/Capstone-Project-Software-Team/davisLab_Software/uploaded_pdbs/6Y4O_sidechain_fixed_fixed.pdb")
# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the protein strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=['A4', 'A145'])
protein.flexibility['A47'].setLibraryRotamers(osprey.WILD_TYPE, 'ILE', 'LYS', 'SER').addWildTypeRotamers().setContinuous()

protein.flexibility['A51'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['B3615', 'B3638'])
ligand.flexibility['B3637'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

        
        
# make the conf space for the protein
proteinConfSpace = osprey.ConfSpace(protein)

# make the conf space for the ligand
ligandConfSpace = osprey.ConfSpace(ligand)

# make the conf space for the protein+ligand complex
complexConfSpace = osprey.ConfSpace([protein, ligand])

# how should we compute energies of molecules?
# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
parallelism = osprey.Parallelism(cpuCores=4,  gpus=0, streamsPerGpu=0)
minimizingEcalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism, isMinimizing=True)

# BBK* needs a rigid energy calculator too, for multi-sequence bounds on K*
rigidEcalc = osprey.SharedEnergyCalculator(minimizingEcalc, isMinimizing=False)


# configure BBK*
bbkstar = osprey.BBKStar(
    proteinConfSpace,
    ligandConfSpace,
    complexConfSpace,
    numBestSequences=4,
    writeSequencesToConsole=True,
    writeSequencesToFile='bbkstar_results_out_A47E.tsv',
    epsilon=0.99,
)

# configure BBK* inputs for each conf space
for info in bbkstar.confSpaceInfos():

	# how should we define energies of conformations?
	eref = osprey.ReferenceEnergies(info.confSpace, minimizingEcalc)
	info.confEcalcMinimized = osprey.ConfEnergyCalculator(info.confSpace, minimizingEcalc, referenceEnergies=eref)

	# compute the energy matrix
	ematMinimized = osprey.EnergyMatrix(info.confEcalcMinimized, cacheFile='emat.%s.dat' % info.id)

	# how should confs be ordered and searched?
	# (since we're in a loop, need capture variables above by using defaulted arguments)
	def makeAStarMinimized(rcs, emat=ematMinimized):
		return osprey.AStarTraditional(emat, rcs, showProgress=False)
	info.confSearchFactoryMinimized = osprey.BBKStar.ConfSearchFactory(makeAStarMinimized)

	# BBK* needs rigid energies too
	confEcalcRigid = osprey.ConfEnergyCalculatorCopy(info.confEcalcMinimized, rigidEcalc)
	ematRigid = osprey.EnergyMatrix(confEcalcRigid, cacheFile='emat.%s.rigid.dat' % info.id)
	def makeAStarRigid(rcs, emat=ematRigid):
		return osprey.AStarTraditional(emat, rcs, showProgress=False)
	info.confSearchFactoryRigid = osprey.BBKStar.ConfSearchFactory(makeAStarRigid)

	# how should we score each sequence?
	# (since we're in a loop, need capture variables above by using defaulted arguments)
	def makePfunc(rcs, confEcalc=info.confEcalcMinimized, emat=ematMinimized):
		return osprey.PartitionFunction(
			confEcalc,
			osprey.AStarTraditional(emat, rcs, showProgress=False),
			osprey.AStarTraditional(emat, rcs, showProgress=False),
			rcs
		)
	info.pfuncFactory = osprey.KStar.PfuncFactory(makePfunc)

# run BBK*
scoredSequences = bbkstar.run(minimizingEcalc.tasks)

# make a sequence analyzer to look at the results
analyzer = osprey.SequenceAnalyzer(bbkstar)

# use results
for scoredSequence in scoredSequences:
	print("result:")
	print("	sequence: %s" % scoredSequence.sequence)
	print("	K* score: %s" % scoredSequence.score)

	# write the sequence ensemble, with up to 10 of the lowest-energy conformations
	numConfs = 10
	analysis = analyzer.analyze(scoredSequence.sequence, numConfs)
	print(analysis)
	analysis.writePdb(
		'seq.%s.pdb' % scoredSequence.sequence,
		'Top %d conformations for sequence %s' % (numConfs, scoredSequence.sequence)
	)
