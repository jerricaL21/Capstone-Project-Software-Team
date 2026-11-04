import os
import re
import json
import pandas as pd

def extract_scores_from_stdout(stdout_text):
    """
    Parses OSPREY stdout text using regex to extract sequence mutations and K* scores.
    Reconstructs WT_Full as 'A38=SER' so logic6.py's parse_residue_site can parse it.
    """
    if not stdout_text:
        return None

    # Matches lines like: sequence  1/  4   A38=ASP   K*(log10): 23.239330 ...
    pattern = re.compile(r"sequence\s+\d+/\s+\d+\s+([A-Za-z0-9_=]+)\s+K\*\([^\)]*\):\s*([\d\.\-]+)")
    matches = pattern.findall(stdout_text)
    if not matches:
        return None

    # Identify wild-type sequence entry (first entry in OSPREY output)
    wt_seq_str, wt_score_str = matches[0]
    wt_score = float(wt_score_str)

    # Parse wild-type residue and amino acid (e.g. "A38" and "ser")
    wt_parts = wt_seq_str.split('=')
    wt_site = wt_parts[0]  # e.g., "A38"
    wt_aa = wt_parts[1].upper()  # e.g., "SER"
    
    # Standardize label to match logic6's expected format (e.g. "A38=SER")
    wt_full_label = f"{wt_site}={wt_aa}"

    parsed_entries = []
    for seq_str, score_str in matches:
        parts = seq_str.split('=')
        if len(parts) == 2:
            mut_aa = parts[1].upper()
            score = float(score_str)
            delta = score - wt_score
            
            parsed_entries.append({
                'WT_Full': wt_full_label,
                'Mutation': mut_aa,
                'Score': round(score, 3),
                'Delta': round(delta, 3)
            })

    return parsed_entries

def read_data_from_json(json_input):
    """Reads execution results list or loads JSON from path."""
    if isinstance(json_input, str) and os.path.exists(json_input):
        with open(json_input, 'r', encoding='utf-8') as f:
            data = json.load(f)
    elif isinstance(json_input, list):
        data = json_input
    else:
        return []

    all_results = []
    for item in data:
        stdout = item.get('stdout', '') if isinstance(item, dict) else ''
        scores = extract_scores_from_stdout(stdout)
        if scores:
            all_results.extend(scores)

    return all_results

def analyze_results(json_input):
    """
    Main entry point for gui.py.
    Processes stdout data and returns (pivot_delta, pivot_annot, summary_df, error_msg).
    """
    if isinstance(json_input, str) and not os.path.exists(json_input):
        all_rows = extract_scores_from_stdout(json_input)
    else:
        all_rows = read_data_from_json(json_input)

    if not all_rows:
        return None, None, None, "No valid K* scores found in stdout data."

    # Create Master DataFrame
    master_df = pd.DataFrame(all_rows)

    # DataFrame 1: Heatmap pivot matrix (0.000 delta = WT reference)
    pivot_delta = master_df.pivot_table(
        index=['WT_Full'], 
        columns='Mutation', 
        values='Delta',
        aggfunc='first'
    ).sort_index()

    # DataFrame 1b: Annotation matrix
    presence = master_df.pivot_table(
        index=['WT_Full'],
        columns='Mutation',
        values='Delta',
        aggfunc='count'
    ).reindex(index=pivot_delta.index, columns=pivot_delta.columns)

    pivot_annot = presence.apply(
        lambda col: col.apply(lambda n: col.name if pd.notna(n) and n > 0 else '')
    )

    # DataFrame 2: Summary table
    summary_df = master_df.copy()
    summary_df['Label'] = summary_df['Mutation'] + " (" + summary_df['Score'].astype(str) + ")"
    mutation_summary = summary_df.groupby(['WT_Full'])['Label'].apply(lambda x: ', '.join(x)).reset_index()
    mutation_summary.rename(columns={'WT_Full': 'Residue Site', 'Label': 'Mutations & Log10 Scores'}, inplace=True)

    return pivot_delta, pivot_annot, mutation_summary, None