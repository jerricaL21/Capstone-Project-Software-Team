# this script reads jsons and creates python scripts for Osprey in individual folder

import os
import json

def create_folders_and_files(json_file, pdb_path, epsilon=0.01, cpu_cores=16):
    # Extract the name of the JSON file (without the extension) to be used as the folder name
    file_name = os.path.splitext(os.path.basename(json_file))[0]
    output_dir = os.path.join(os.path.dirname(json_file), file_name)
    os.makedirs(output_dir, exist_ok=True)

    # Read the JSON data from the file
    with open(json_file, 'r') as f:
        data = json.load(f)

    # Check if the key "Complex_Size" exists in the JSON data
    if "Complex_Size" in data:
        # Extract the variables from the "Complex_Size" dictionary
        variables = data["Complex_Size"]
        variable1 = variables.get("mutant_beginning", "KEY1")
        variable2 = variables.get("mutant_ending", "KEY2")
        variable3 = variables.get("ligand_beginning", "KEY3")
        variable4 = variables.get("ligand_ending", "KEY4")
        
    else:
        # If "Complex_Size" key is not present, set default values for variables
        variable1 = "KEY1"
        variable2 = "KEY2"
        variable3 = "KEY3"
        variable4 = "KEY4"
        
    generated_files = []# Making an empty

    # Loop through the JSON data and create folders and .py files
    for key, subdata in data.items():
        
        # Skip processing if the key is "Complex_Size"
        if key == "Complex_Size":
            continue

        # Remove the last character from the key
        truncated_key = key[:-1]
    
        # Create a folder with the name of the key (if it doesn't exist)
        subfolder = os.path.join(output_dir, key)
        os.makedirs(subfolder, exist_ok=True)

        # Get the subkeys for object "A" and "B"
        subkeys_from_A = list(subdata["A"].keys())
        subkeys_from_B = list(subdata["B"].keys())            
            
        # Create the content for the .py file
        py_content = f'''

import osprey
osprey.start(heapSizeMiB=100000, garbageSizeMiB=8192)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('{pdb_path}')
# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the protein strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=['{variable1}', '{variable2}'])
protein.flexibility['{truncated_key}'].setLibraryRotamers(osprey.WILD_TYPE, 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'GLU', 'TYR', 'ASP', 'ASN', 'GLN', 'ARG', 'LYS', 'SER', 'THR', 'CYS', 'HIS').addWildTypeRotamers().setContinuous()

# define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['{variable3}', '{variable4}'])
'''

        # Add subkeys from object "A" to the .py content
        for subkey_A in subkeys_from_A:
            py_content += f"protein.flexibility['{subkey_A}'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()\n"

        # Add subkeys from object "B" to the .py content
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
parallelism = osprey.Parallelism(cpuCores={cpu_cores}, gpus=1, streamsPerGpu=82)
minimizingEcalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism, isMinimizing=True)

# BBK* needs a rigid energy calculator too, for multi-sequence bounds on K*
rigidEcalc = osprey.SharedEnergyCalculator(minimizingEcalc, isMinimizing=False)


# configure BBK*
bbkstar = osprey.BBKStar(
	proteinConfSpace,
	ligandConfSpace,
	complexConfSpace,
	numBestSequences=1, # more sequenses - more computation time
	epsilon=0.01, # you proabably want something more precise in your real designs
	writeSequencesToConsole=True,
	writeSequencesToFile='bbkstar_results_{file_name}_{key}.tsv'
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
        with open(py_file, 'w') as f:
            f.write(py_content)
        generated_files.append(py_file)
    return generated_files

