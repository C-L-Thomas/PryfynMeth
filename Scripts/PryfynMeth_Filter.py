import os
import argparse

# Step 1: Filter based on total coverage threshold
def filter_file(input_path, output_path, threshold):
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        header = infile.readline()
        outfile.write(header)

        for line in infile:
            fields = line.strip().split('\t')
            if not fields or len(fields) < 8:
                continue

            try:
                total = int(fields[5])
            except ValueError:
                total = 0

            if total < threshold:
                fields[3] = fields[4] = fields[5] = '0'
                fields[6] = fields[7] = '1.0'

            outfile.write('\t'.join(fields) + '\n')

# Step 2: Filter to only include significant rows (FDR < 0.05)
def methylation_significant_only(filtered_dir, output_dir_m):
    os.makedirs(output_dir_m, exist_ok=True)

    for filename in os.listdir(filtered_dir):
        if filename.startswith('.'):
            continue
        input_path = os.path.join(filtered_dir, filename)
        output_path = os.path.join(output_dir_m, filename)

        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            header = infile.readline()
            outfile.write(header)

            for line in infile:
                fields = line.strip().split('\t')
                if not fields or len(fields) < 8:
                    continue

                try:
                    fdr = float(fields[7])
                except ValueError:
                    continue

                if fdr < 0.05:
                    outfile.write('\t'.join(fields) + '\n')

# Helper: Extract chr:pos:strand as unique key
def build_significant_shared_set(methylated_dir):
    shared_keys = set()
    for filename in os.listdir(methylated_dir):
        if filename.startswith('.'):
            continue
        path = os.path.join(methylated_dir, filename)
        with open(path, 'r') as infile:
            _ = infile.readline()
            for line in infile:
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    key = f"{fields[0]}:{fields[1]}:{fields[2]}"
                    shared_keys.add(key)
    return shared_keys

# Helper: Build key→row map from a file
def build_key_row_map(file_path):
    key_to_row = {}
    with open(file_path, 'r') as f:
        header = f.readline()
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 8:
                key = f"{fields[0]}:{fields[1]}:{fields[2]}"
                key_to_row[key] = fields
    return header, key_to_row

# Step 3: Generate shared methylation output
def write_shared_table(shared_keys, methylated_dir, shared_dir, input_dir, threshold, revert=False):
    os.makedirs(shared_dir, exist_ok=True)

    # If revert is active, build original row maps
    original_data = {}
    if revert:
        for filename in os.listdir(input_dir):
            if filename.startswith('.'):
                continue
            input_path = os.path.join(input_dir, filename)
            _, key_to_row = build_key_row_map(input_path)
            original_data[filename] = key_to_row

    for filename in os.listdir(methylated_dir):
        if filename.startswith('.'):
            continue

        methylated_path = os.path.join(methylated_dir, filename)
        shared_path = os.path.join(shared_dir, filename)

        header, methylated_rows = build_key_row_map(methylated_path)
        original_rows = original_data.get(filename, {}) if revert else {}

        with open(shared_path, 'w') as out:
            out.write(header)
            for key in sorted(shared_keys, key=lambda k: (k.split(':')[0], int(k.split(':')[1]), k.split(':')[2])):
                if key in methylated_rows:
                    row = methylated_rows[key]
                elif revert and key in original_rows:
                    original = original_rows[key]
                    try:
                        original_total = int(original[5])
                    except ValueError:
                        original_total = 0

                    if original_total >= threshold:
                        row = original
                    else:
                        row = [original[0], original[1], original[2], '0', '0', '0', '1.0', '1.0']
                else:
                    chrom, pos, strand = key.split(':')
                    row = [chrom, pos, strand, '0', '0', '0', '1.0', '1.0']
                out.write('\t'.join(row) + '\n')

# Main control
def main():
    parser = argparse.ArgumentParser(description="Filter and merge methylation call files.")
    parser.add_argument("-i", "--input", required=True, help="Input folder path")
    parser.add_argument("-f", "--filtered", required=True, help="Filtered output folder")
    parser.add_argument("-m", "--methylated", required=True, help="Methylated significant output folder")
    parser.add_argument("-s", "--shared", required=True, help="Shared methylation output folder")
    parser.add_argument("-threshold", type=int, required=True, help="Coverage threshold")
    parser.add_argument("-revert", action="store_true", help="Restore original values for non-significant sites in shared output, if above threshold")

    args = parser.parse_args()

    os.makedirs(args.filtered, exist_ok=True)

    print("[1/3] Filtering by total coverage...")
    for filename in os.listdir(args.input):
        if filename.startswith('.'):
            continue
        input_file = os.path.join(args.input, filename)
        filtered_file = os.path.join(args.filtered, filename)
        filter_file(input_file, filtered_file, args.threshold)

    print("[2/3] Keeping rows with FDR < 0.05...")
    methylation_significant_only(args.filtered, args.methylated)

    print("[3/3] Building shared significant site table...")
    shared_keys = build_significant_shared_set(args.methylated)
    write_shared_table(shared_keys, args.methylated, args.shared, args.input, args.threshold, revert=args.revert)

    print("✅ Done!")

if __name__ == "__main__":
    main()
