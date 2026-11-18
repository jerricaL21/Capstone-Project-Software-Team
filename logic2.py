#this script sort keys by the distances for each object (A and B) and add variables to json
import os
import json

def add_variables_to_json(json_file, mutant_beginning, mutant_ending, ligand_beginning, ligand_ending):
    # Read the JSON data from the file
    with open(json_file, 'r') as f:
        data = json.load(f)
        
     # Build the dictionary internally
    data["Complex_Size"] = {
        "mutant_beginning": mutant_beginning,
        "mutant_ending":    mutant_ending,
        "ligand_beginning": ligand_beginning,
        "ligand_ending":    ligand_ending,
    }

    # Sort the residues inside each chain object (any chain ID, not just A and B)
    # from closest to farthest. Only lists of numeric distances are sorted, so
    # "Complex_Size" and anything non-numeric is left untouched.
    for key in data:
        if key == "Complex_Size" or not isinstance(data[key], dict):
            continue
        for chain_id, contacts in data[key].items():
            if (isinstance(contacts, dict) and contacts
                    and all(isinstance(v, (int, float)) for v in contacts.values())):
                data[key][chain_id] = dict(sorted(contacts.items(), key=lambda item: item[1]))
        
    # Save the modified JSON with the same name
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=4)
