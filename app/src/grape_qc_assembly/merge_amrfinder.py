#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd

def get_sort_key(scope, type_str):
    """
    Returns a sorting key for Scope and Type combinations.
    - Scope priority: core > plus > others alphabetically (case-insensitive)
    - Type priority: AMR > virulence > STRESS > others alphabetically (case-insensitive)
    """
    s = str(scope).strip().lower()
    t = str(type_str).strip().lower()
    
    # Scope priority: core > plus > others alphabetically
    if s == 'core':
        scope_pri = 0
    elif s == 'plus':
        scope_pri = 1
    else:
        scope_pri = 2
        
    # Type priority: amr > virulence > stress > others alphabetically
    if t == 'amr':
        type_pri = 0
    elif t == 'virulence':
        type_pri = 1
    elif t == 'stress':
        type_pri = 2
    else:
        type_pri = 3
        
    return (scope_pri, s, type_pri, t)

def get_group_name(scope, type_str):
    """
    Formats the group header name, e.g. "core_AMR", "plus_virulence"
    """
    s = str(scope).strip().lower()
    t = str(type_str).strip().lower()
    
    if t == 'amr':
        t_norm = 'AMR'
    elif t == 'virulence':
        t_norm = 'virulence'
    elif t == 'stress':
        t_norm = 'stress'
    else:
        t_norm = t
        
    return f"{s}_{t_norm}"

def merge_amrfinder_files(input_dir, list_path, output_dir):
    """
    Core logic to merge amrfinder-plus output TSV files.
    Raises standard exceptions for error conditions:
    - FileNotFoundError: If the strain list file is missing.
    - ValueError: If the strain list is empty or no valid TSV files could be parsed.
    - OSError: For directory creation or file writing errors.
    """
    # 1. Read strain list
    if not os.path.exists(list_path):
        raise FileNotFoundError(f"Strain list file '{list_path}' does not exist.")
        
    strains = []
    try:
        with open(list_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line:
                    strains.append(line)
    except (OSError, UnicodeDecodeError) as e:
        raise OSError(f"Failed to read strain list file '{list_path}': {e}")
                
    if not strains:
        raise ValueError(f"No strain IDs found in the list file '{list_path}'.")
        
    # Ensure the output directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as e:
        raise OSError(f"Failed to create output directory '{output_dir}': {e}")
    
    dfs = []
    
    # 2. Load and parse TSV for each strain
    for strain in strains:
        tsv_path = None
        candidates = [
            f"{strain}_plus_out.tsv",
            f"{strain}_out.tsv"
        ]
        for c in candidates:
            p = os.path.join(input_dir, c)
            if os.path.exists(p):
                tsv_path = p
                break
                
        if not tsv_path:
            raise FileNotFoundError(f"TSV file for strain '{strain}' not found in '{input_dir}'.")
            
        try:
            df = pd.read_csv(tsv_path, sep='\t')
            # Normalize column names to lowercase/trimmed to handle variation
            df.columns = [c.strip().lower() for c in df.columns]
            
            required = ['element symbol', 'scope', 'type']
            if not all(col in df.columns for col in required):
                missing = [col for col in required if col not in df.columns]
                raise KeyError(f"Missing required columns: {missing}")
                
            # Filter rows with valid element symbols
            df = df[required].dropna(subset=['element symbol'])
            df = df.rename(columns={'element symbol': 'gene'})
            df['strain'] = strain
            dfs.append(df)
        except pd.errors.EmptyDataError as e:
            raise ValueError(f"TSV file '{tsv_path}' is empty: {e}")
        except pd.errors.ParserError as e:
            raise ValueError(f"TSV file '{tsv_path}' is malformed: {e}")
        except KeyError as e:
            raise KeyError(f"Column error in '{tsv_path}': {e}")
        except (OSError, UnicodeDecodeError) as e:
            raise OSError(f"Could not read TSV file '{tsv_path}': {e}")
            
    # Combine all data
    combined_df = pd.concat(dfs, ignore_index=True)
    # Deduplicate genes per strain (presence-absence is binary)
    combined_df = combined_df.drop_duplicates(subset=['strain', 'gene'])
    
    # Deduplicate genes to find their unique group mapping
    gene_info = combined_df[['gene', 'scope', 'type']].drop_duplicates(subset=['gene']).copy()
    gene_info['group'] = gene_info.apply(
        lambda row: get_group_name(row['scope'], row['type']), axis=1
    )
    
    # 3. Sort groups
    unique_groups = gene_info[['group', 'scope', 'type']].drop_duplicates().copy()
    unique_groups['sort_key'] = unique_groups.apply(
        lambda row: get_sort_key(row['scope'], row['type']), axis=1
    )
    sorted_groups = unique_groups.sort_values(by='sort_key')['group'].tolist()
    
    # 4. Construct final column layout
    columns = ['Strain']
    separator_cols = set(sorted_groups)
    
    for grp in sorted_groups:
        columns.append(grp)
        # get genes belonging to this group, sorted alphabetically
        grp_genes = sorted(gene_info[gene_info['group'] == grp]['gene'].unique())
        columns.extend(grp_genes)
        
    # 5. Populate rows for the final matrix
    matrix_rows = []
    for strain in strains:
        row = {'Strain': strain}
        present_genes = set(combined_df[combined_df['strain'] == strain]['gene'])
        
        for col in columns[1:]:
            if col in separator_cols:
                row[col] = '|'
            else:
                row[col] = 1 if col in present_genes else 0
        matrix_rows.append(row)
        
    # Create the final summary DataFrame
    summary_df = pd.DataFrame(matrix_rows, columns=columns)
    
    # 6. Determine output file path and save
    list_basename = os.path.basename(list_path)
    list_name, _ = os.path.splitext(list_basename)
    output_filename = f"{list_name}_summary.tsv"
    output_path = os.path.join(output_dir, output_filename)
    
    try:
        summary_df.to_csv(output_path, sep='\t', index=False)
        print(f"Successfully generated summary TSV at: {output_path}")
    except OSError as e:
        raise OSError(f"Failed to write output file '{output_path}': {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Merge amrfinder-plus output TSV files using pandas."
    )
    parser.add_argument(
        "input_dir", 
        help="Path to the folder containing amrfinder TSV files."
    )
    parser.add_argument(
        "list_path", 
        help="Path to the strain list text file (e.g., list_amrfinder.txt)."
    )
    parser.add_argument(
        "output_dir", 
        help="Path to the folder where the merged summary TSV will be saved."
    )
    
    args = parser.parse_args()
    
    try:
        merge_amrfinder_files(args.input_dir, args.list_path, args.output_dir)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except OSError as e:
        print(f"OS Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
