import os
import sys
import argparse
import scipy.stats as stats

def perform_binomial_test(successes, total, lambda_value):
    return stats.binomtest(successes, total, lambda_value, alternative='two-sided').pvalue

def fdr_correction(p_values):
    sorted_pvals = sorted((p, i) for i, p in enumerate(p_values))
    m = len(p_values)
    corrected_pvals = [0] * m
    for rank, (p, i) in enumerate(sorted_pvals, start=1):
        corrected_pvals[i] = min(p * m / rank, 1.0)
    return corrected_pvals

def read_metadata(file_path):
    metadata = []
    with open(file_path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            metadata.append(dict(zip(header, fields)))
    return metadata

def read_input_illu(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            data.append(fields)
    return data

def read_input_nano(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chr, pos, strand, C, T, total = fields[0], fields[1], fields[5], int(fields[11]), int(fields[12]), int(fields[9])
            data.append((chr, pos, strand, C, T, total))
    return data

def run_analysis(metadata_file, platform, input_dir, output_dir):
    if platform not in ['illu', 'nano']:
        raise ValueError("Unsupported platform specified. Use 'illu' or 'nano'.")

    if not os.path.isfile(metadata_file):
        raise FileNotFoundError(f"Metadata file '{metadata_file}' does not exist.")
    if not os.path.isdir(input_dir):
        raise NotADirectoryError(f"Input directory '{input_dir}' does not exist.")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    metadata = read_metadata(metadata_file)

    for row in metadata:
        sample_name = row['sampleName']
        input_file = os.path.join(input_dir, row['inputFile'])
        lambda_value = float(row['lambda'])

        if not os.path.isfile(input_file):
            print(f"Input file '{input_file}' for sample '{sample_name}' not found. Skipping.")
            continue

        output_file = os.path.join(output_dir, f"{sample_name}_binomial_results.txt")

        try:
            if platform == 'illu':
                input_data = read_input_illu(input_file)
            else:
                input_data = read_input_nano(input_file)
        except Exception as e:
            print(f"Error reading input file '{input_file}': {e}")
            continue

        results = []
        p_values = []

        for fields in input_data:
            if platform == 'illu':
                chr, pos, strand, C, T, _, _ = fields
                C = int(C)
                T = int(T)
                total = C + T
            else:
                chr, pos, strand, C, T, total = fields

            p_value = 1.0 if total < 1 else perform_binomial_test(C, total, lambda_value)
            p_values.append(p_value)
            results.append((chr, pos, strand, C, T, total, p_value))

        corrected_pvals = fdr_correction(p_values)

        with open(output_file, 'w') as out_f:
            out_f.write("chr\tpos\tstrand\tC\tT\ttotal\tp-value\tfdr\n")
            for result, fdr_pval in zip(results, corrected_pvals):
                out_f.write("\t".join(map(str, result)) + f"\t{fdr_pval}\n")

    print("Processing complete. Results are saved in the specified output directory.")

def parse_args():
    parser = argparse.ArgumentParser(description="Run binomial methylation test.")
    parser.add_argument("-meta", required=True, help="Path to metadata file")
    parser.add_argument("-platform", required=True, choices=["illu", "nano"], help="Platform (illu or nano)")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
    return parser.parse_args()

def main():
    args = parse_args()
    run_analysis(
        metadata_file=args.meta,
        platform=args.platform,
        input_dir=args.input_dir,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main()
