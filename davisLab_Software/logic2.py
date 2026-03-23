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

    # Sort the keys inside each object (A and B) based on their values
    for key in data:
        if isinstance(data[key], dict):
            if "A" in data[key]:
                data[key]["A"] = {k: v for k, v in sorted(data[key]["A"].items(), key=lambda item: item[1])}
            if "B" in data[key]:
                data[key]["B"] = {k: v for k, v in sorted(data[key]["B"].items(), key=lambda item: item[1])}
        
    # Save the modified JSON with the same name
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=4)


