# this script reads jsons and creates python scripts for Osprey in individual folder

import os
import json
import sys
import subprocess

def create_folders_and_files(json_file, pdb_path, epsilon=0.01, cpu_cores=16, gpus = 0, streams_per_gpu = 0, heap_size = 100000, garbage_size = 8192, target_residue=None, amino_acids=None):
    # Extract the name of the JSON file (without the extension) to be used as the folder name
    pdb_path  = pdb_path.replace("\\", "/")
    json_file = json_file.replace("\\", "/")
    file_name = os.path.splitext(os.path.basename(json_file))[0]
    output_dir = os.path.join(os.path.dirname(json_file), file_name)
    os.makedirs(output_dir, exist_ok=True)

    # Read the JSON data from the file
    with open(json_file, 'r') as f:
        data = json.load(f)

    # Validate Complex_Size exists and is populated
    if "Complex_Size" not in data:
        raise ValueError("JSON is missing 'Complex_Size' block. Please complete Step 2 first.")

    variables = data["Complex_Size"]
    variable1 = variables.get("mutant_beginning") or ""
    variable2 = variables.get("mutant_ending") or ""
    variable3 = variables.get("ligand_beginning") or ""
    variable4 = variables.get("ligand_ending") or ""

    if not all([variable1, variable2, variable3, variable4]):
        raise ValueError(f"Complex_Size has empty fields: {variables}. Please re-run Step 2 with valid residue ranges.")


    #-----Requires at least one amino acid to be selected. If not, OSPREY does't run---#
    if not amino_acids:
        raise ValueError("Please select at least one amino acid before generating OSPREY scripts.")    
    
    generated_files = []# Making an empty

    # Loop through the JSON data and create folders and .py files
    for key, subdata in data.items():
        
        # Skip processing if the key is "Complex_Size"
        if key == "Complex_Size":
            continue
        # Only process the selected residue
        if target_residue and key != target_residue:
            continue
        # Remove the last character from the key
        # Format: 'A 7' from key 'A7E'
        chain = key[0]
        number = key[1:-1]
        truncated_key = f"{chain}{number}"
    
        # Create a folder with the name of the key (if it doesn't exist)
        subfolder = os.path.join(output_dir, key)
        os.makedirs(subfolder, exist_ok=True)

        # Get the subkeys for object "A" and "B"
        subkeys_from_A = list(subdata["A"].keys())
        subkeys_from_B = list(subdata["B"].keys())            
            
        # Create the content for the .py file

        # Build amino acid string from user selection
        if amino_acids:
            aa_string = ", ".join(f"'{aa}'" for aa in amino_acids)
            # The count is the number of mutations + 1 for the WILD_TYPE
            num_seqs = len(amino_acids) + 1
        else:
            aa_string = "'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'GLU', 'TYR', 'ASP', 'ARG', 'ASN', 'CYS', 'GLN', 'GLY', 'HIS', 'LYS', 'SER', 'THR'"
            # There are 18 amino acids above + 1 for the WILD_TYPE
            num_seqs = 11
        py_content = f'''

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

    candidates.extend(glob.glob(r"C:\\Program Files\\Java\\jdk-*"))
    candidates.extend(glob.glob(r"C:\\Program Files\\Eclipse Adoptium\\jdk-*"))

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
        heapSizeMiB={heap_size},
        enableAssertions=False,
        stackSizeMiB=16,
        garbageSizeMiB={garbage_size},
        allowRemoteManagement=False,
        attachJvmDebugger=False,
    ))
    if osprey._print_preamble:
        print("Using up to {heap_size} MiB heap memory: {garbage_size} MiB for garbage, %d MiB for storage" % ({heap_size} - {garbage_size}))
else:
    # No system JDK found on this machine - fall back to OSPREY's bundled JRE
    osprey.start(heapSizeMiB={heap_size}, garbageSizeMiB={garbage_size})

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb("{pdb_path}")
# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the protein strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=['{variable1}', '{variable2}'])
protein.flexibility['{truncated_key}'].setLibraryRotamers(osprey.WILD_TYPE, {aa_string}).addWildTypeRotamers().setContinuous()

'''

        # Add subkeys from object "A" to the .py content
        for subkey_A in subkeys_from_A:
            py_content += f"protein.flexibility['{subkey_A}'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()\n"

        py_content += f"\n# define the ligand strand\nligand = osprey.Strand(mol, templateLib=templateLib, residues=['{variable3}', '{variable4}'])\n"


        for subkey_B in subkeys_from_B:
            py_content += f"ligand.flexibility['{subkey_B}'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()\n"
        py_content += f'''
        
        
# make the conf space for the protein
proteinConfSpace = osprey.ConfSpace(protein)

# make the conf space for the ligand
ligandConfSpace = osprey.ConfSpace(ligand)

# make the conf space for the protein+ligand complex
complexConfSpace = osprey.ConfSpace([protein, ligand])

# how should we compute energies of molecules?
# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
parallelism = osprey.Parallelism(cpuCores={cpu_cores},  gpus={gpus}, streamsPerGpu={streams_per_gpu})
minimizingEcalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism, isMinimizing=True)

# BBK* needs a rigid energy calculator too, for multi-sequence bounds on K*
rigidEcalc = osprey.SharedEnergyCalculator(minimizingEcalc, isMinimizing=False)


# configure BBK*
bbkstar = osprey.BBKStar(
    proteinConfSpace,
    ligandConfSpace,
    complexConfSpace,
    numBestSequences={num_seqs},
    writeSequencesToConsole=True,
    writeSequencesToFile='bbkstar_results_{file_name}_{key}.tsv',
    epsilon={epsilon},
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
	print("\tsequence: %s" % scoredSequence.sequence)
	print("\tK* score: %s" % scoredSequence.score)

	# write the sequence ensemble, with up to 10 of the lowest-energy conformations
	numConfs = 10
	analysis = analyzer.analyze(scoredSequence.sequence, numConfs)
	print(analysis)
	analysis.writePdb(
		'seq.%s.pdb' % scoredSequence.sequence,
		'Top %d conformations for sequence %s' % (numConfs, scoredSequence.sequence)
	)
'''
        # Write the .py file inside the key's folder
        py_file = os.path.join(subfolder, f'bbkstar_{file_name}_{key}.py')
        with open(py_file, 'w', encoding='utf-8') as f:
            f.write(py_content)
        generated_files.append(py_file)
    return generated_files

def run_osprey_scripts(generated_files):
    """Run all generated bbkstar*.py scripts using the current Python interpreter"""
    results = []
    for path in generated_files:
        folder = os.path.dirname(path)
        try:
            result = subprocess.run(
                [sys.executable, path],  # ← uses whatever python is running the app
                cwd=folder,
                capture_output=True,
                text=True
            )
            results.append({
                "script": path,
                "returncode": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr
            })
        except Exception as e:
            results.append({"script": path, "error": str(e)})
    return results