import os
import sys
import scipy.stats as stats

# Function to perform one-sided binomial test (testing if methylation > background)
def perform_binomial_test(successes, total, lambda_value):
    # One-sided alternative: test if observed methylation rate is greater than lambda background rate
    return stats.binomtest(successes, total, lambda_value, alternative='greater').pvalue

# Function to perform FDR correction using Benjamini-Hochberg
def fdr_correction(p_values):
    sorted_pvals = sorted((p, i) for i, p in enumerate(p_values))
    m = len(p_values)
    corrected_pvals = [0] * m
    for rank, (p, i) in enumerate(sorted_pvals, start=1):
        corrected_pvals[i] = min(p * m / rank, 1.0)
    return corrected_pvals

# Function to read metadata
def read_metadata(file_path):
    metadata = []
    with open(file_path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            metadata.append(dict(zip(header, fields)))
    return metadata

# Function to read input file for the 'illu' platform
def read_input_illu(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            data.append(fields)
    return data

# ✅ Corrected function to read 'nano' platform input
def read_input_nano(file_path):
    data = []
    with open(file_path, 'r') as f:
        for i, line in enumerate(f, 1):
            fields = line.strip().split('\t')
            if len(fields) < 13:
                print(f"[WARN] Skipping line {i} in '{file_path}': only {len(fields)} fields (expected ≥13)")
                continue
            try:
                chr = fields[0]
                pos = fields[1]
                strand = fields[5]  # ✅ use column 6 for strand
                total = int(fields[4])     # Column 5 (0-based index 4)
                C = int(fields[11])        # Column 12 (0-based index 11)
                T = int(fields[12])        # Column 13 (0-based index 12)
                data.append((chr, pos, strand, C, T, total))
            except Exception as e:
                print(f"[ERROR] Failed to parse line {i} in '{file_path}': {e}")
                continue
    print(f"[INFO] Read {len(data)} valid rows from '{file_path}'")
    return data

# Main script
if len(sys.argv) < 9:
    print("Error: Insufficient arguments provided.")
    print("Usage: python3 PryfMeth_Binomial.py -meta <metadata_file> -platform <platform> -i <input_dir> -o <output_dir>")
    sys.exit(1)

# Parse command-line arguments
args = sys.argv[1:]
metadata_file = None
platform = None
input_dir = None
output_dir = None

for i in range(0, len(args), 2):
    if args[i] == '-meta':
        metadata_file = args[i + 1]
    elif args[i] == '-platform':
        platform = args[i + 1]
    elif args[i] == '-i':
        input_dir = args[i + 1]
    elif args[i] == '-o':
        output_dir = args[i + 1]

if not metadata_file or not platform or not input_dir or not output_dir:
    print("Error: Missing required arguments.")
    print("Usage: python3 PryfMeth_Binomial.py -meta <metadata_file> -platform <platform> -i <input_dir> -o <output_dir>")
    sys.exit(1)

if platform not in ['illu', 'nano']:
    print("Error: Unsupported platform specified. Use 'illu' or 'nano'.")
    sys.exit(1)

if not os.path.isfile(metadata_file):
    print(f"Error: Metadata file '{metadata_file}' does not exist.")
    sys.exit(1)

if not os.path.isdir(input_dir):
    print(f"Error: Input directory '{input_dir}' does not exist.")
    sys.exit(1)

if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

# Read the metadata
metadata = read_metadata(metadata_file)

# Process each sample
for row in metadata:
    sample_name = row['sampleName']
    input_file = os.path.join(input_dir, row['inputFile'])
    lambda_value = float(row['lambda'])

    if not os.path.isfile(input_file):
        print(f"[WARN] Input file '{input_file}' for sample '{sample_name}' not found. Skipping.")
        continue

    output_file = os.path.join(output_dir, f"{sample_name}_binomial_results.txt")

    try:
        if platform == 'illu':
            input_data = read_input_illu(input_file)
        elif platform == 'nano':
            input_data = read_input_nano(input_file)
    except Exception as e:
        print(f"[ERROR] Failed to read input file '{input_file}': {e}")
        continue

    results = []
    p_values = []

    for fields in input_data:
        if platform == 'illu':
            chr, pos, strand, C, T, _, _ = fields
            C = int(C)
            T = int(T)
            total = C + T
        elif platform == 'nano':
            chr, pos, strand, C, T, total = fields

        if total < 1:
            p_value = 1.0
        else:
            p_value = perform_binomial_test(C, total, lambda_value)

        p_values.append(p_value)
        results.append((chr, pos, strand, C, T, total, p_value))

    if not results:
        print(f"[WARN] No valid data rows processed for sample '{sample_name}'. Output will be empty.")
        continue

    corrected_pvals = fdr_correction(p_values)

    with open(output_file, 'w') as out_f:
        out_f.write("chr\tpos\tstrand\tC\tT\ttotal\tp-value\tfdr\n")
        for result, fdr_pval in zip(results, corrected_pvals):
            out_f.write("\t".join(map(str, result)) + f"\t{fdr_pval}\n")

    print(f"[INFO] Processed {len(results)} rows for sample '{sample_name}'. Output saved to '{output_file}'.")

print("✅ Processing complete (one-sided binomial test used).")
