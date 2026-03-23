import osprey
osprey.start()

ffparams = osprey.ForcefieldParams()
mol = osprey.readPdb("uploaded_pdbs/1IWQ_fixed.pdb")
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

protein = osprey.Strand(mol, templateLib=templateLib, residues=['A7', 'A146'])

print("Testing flexibility lookup:")
res = protein.flexibility['A7']
print(f"A7 flexibility: {res}")