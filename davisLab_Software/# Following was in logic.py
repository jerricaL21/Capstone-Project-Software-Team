# Following was in logic.py after the import section, the first function in the file



def fix_pdb_protonation(pdb_path):
    fixed_path = pdb_path.replace(".pdb", "_fixed.pdb")
    
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            res_name = line[17:20].strip()
            # Skip water only
            if res_name in ["HOH", "WAT"]:
                continue
            # Fix histidine protonation
            if res_name in ["HIS", "HISA", "HISD", "HISE", "HISH"]:
                line = line[:17] + "HIE" + line[20:]
        fixed_lines.append(line)
    
    with open(fixed_path, 'w') as f:
        f.writelines(fixed_lines)
    
    return fixed_path