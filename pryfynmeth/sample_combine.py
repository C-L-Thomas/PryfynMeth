import os
import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge methylation files grouped by condition using metadata.")
    parser.add_argument("-i", "--input", required=True, help="Input folder containing methylation files")
    parser.add_argument("-m", "--metadata", required=True, help="CSV file with 'fileName' and 'condition' columns")
    parser.add_argument("-o", "--output", required=True, help="Output folder for merged condition files")
    parser.add_argument("-t", "--type", choices=['nano', 'illu'], required=True, help="Input type: 'nano' or 'illu'")
    return parser.parse_args()

def read_nano_file(filepath, min_cols):
    df = pd.read_csv(filepath, sep="\t", header=None, low_memory=False)
    for i in range(df.shape[1], min_cols):
        df[i] = 0

    for col in [4, 11, 12]:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

    # FIX: Correct strand column is index 5 (column 6 in input)
    df[5] = df[5].astype(str).str.strip()
    df[5] = df[5].where(df[5].isin(['+', '-']), '-')

    # Return with corrected strand column
    return df[[0, 1, 2, 5, 4, 11, 12]]  # chr, start, end, strand, coverage, meth, unmeth

def read_illu_file(filepath):
    df = pd.read_csv(filepath, sep="\t", header=None, low_memory=False)
    for col in [3, 4]:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    return df[[0, 1, 2, 3, 4, 5, 6]]  # chr, pos, strand, meth, unmeth, context, sequence

def merge_nano(file_paths):
    min_cols = 14
    merged = read_nano_file(file_paths[0], min_cols)
    merged.columns = ['chr', 'start', 'end', 'strand', 'coverage', 'meth', 'unmeth']

    for path in file_paths[1:]:
        df = read_nano_file(path, min_cols)
        df.columns = ['chr', 'start', 'end', 'strand', 'coverage_new', 'meth_new', 'unmeth_new']
        merged = pd.merge(merged, df, on=['chr', 'start', 'end', 'strand'], how='outer')

        for col in ['coverage', 'meth', 'unmeth']:
            merged[col] = merged[col].fillna(0) + merged[f"{col}_new"].fillna(0)
            merged.drop(columns=[f"{col}_new"], inplace=True)

    merged['percent'] = (
        merged['meth'] / merged['coverage']
    ).replace([float('inf'), -float('inf')], 0).fillna(0) * 100
    merged['percent'] = merged['percent'].round(2)

    merged['meth_status'] = '-'
    merged['start2'] = merged['start']
    merged['end2'] = merged['end']
    for col in ['col9', 'col10', 'col15', 'col16', 'col17', 'col18', 'col19']:
        merged[col] = 0

    merged['coverage'] = merged['coverage'].round().astype(int)

    final_cols = [
        'chr', 'start', 'end', 'meth_status', 'coverage', 'strand',
        'start2', 'end2', 'col9', 'col10', 'percent', 'meth', 'unmeth',
        'col15', 'col16', 'col17', 'col18', 'col19'
    ]

    merged = merged[final_cols]
    merged.sort_values(by=['chr', 'start', 'end'], inplace=True)
    return merged

def merge_illu(file_paths):
    merged = read_illu_file(file_paths[0])
    merged.columns = ['chr', 'pos', 'strand', 'meth', 'unmeth', 'context', 'sequence']

    for path in file_paths[1:]:
        df = read_illu_file(path)
        df.columns = ['chr', 'pos', 'strand', 'meth_new', 'unmeth_new', 'context_new', 'sequence_new']
        merged = pd.merge(merged, df, on=['chr', 'pos', 'strand'], how='outer')

        for col in ['meth', 'unmeth']:
            merged[col] = merged[col].fillna(0) + merged[f"{col}_new"].fillna(0)
            merged.drop(columns=[f"{col}_new"], inplace=True)

        # Drop new context/sequence — keep original
        merged.drop(columns=['context_new', 'sequence_new'], inplace=True)

    merged = merged[['chr', 'pos', 'strand', 'meth', 'unmeth', 'context', 'sequence']]
    merged.sort_values(by=['chr', 'pos'], inplace=True)
    return merged

def main():
    args = parse_arguments()
    input_dir = args.input
    metadata_file = args.metadata
    output_dir = args.output
    file_type = args.type

    os.makedirs(output_dir, exist_ok=True)

    metadata = pd.read_csv(metadata_file, sep=None, engine='python')
    if 'fileName' not in metadata.columns or 'condition' not in metadata.columns:
        print("Metadata file must contain 'fileName' and 'condition' columns.")
        return

    grouped = metadata.groupby('condition')

    for condition, group in grouped:
        file_list = group['fileName'].tolist()
        file_paths = [os.path.join(input_dir, f) for f in file_list]

        print(f"Processing condition '{condition}' with {len(file_paths)} files...")

        if file_type == 'nano':
            merged_df = merge_nano(file_paths)
        else:
            merged_df = merge_illu(file_paths)

        output_path = os.path.join(output_dir, condition)
        merged_df.to_csv(output_path, sep="\t", header=False, index=False)
        print(f" → Output written to: {output_path}")

if __name__ == "__main__":
    main()
