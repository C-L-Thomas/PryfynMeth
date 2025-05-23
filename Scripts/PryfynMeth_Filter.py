import os
import argparse

# Step 1: Filter based on total coverage and FDR
def filter_file(input_path, output_path, threshold):
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        header = infile.readline()
        outfile.write(header)

        for line in infile:
            fields = line.strip().split('\t')
            if not fields or len(fields) < 8:
                continue

            try:
                coverage = int(fields[5])
                fdr = float(fields[7])  # FDR is in column 8 (index 7)
            except ValueError:
                coverage = 0
                fdr = 1.0

            if coverage < threshold:
                # Zero everything
                fields[3] = fields[4] = fields[5] = '0'
                fields[6] = fields[7] = '1.0'
            elif fdr >= 0.05:
                # High FDR: zero unmethylated count only (column 4, index 3)
                fields[3] = '0'

            outfile.write('\t'.join(fields) + '\n')

# Step 2: Keep only rows with FDR < 0.05
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

# Step 3a: Create combined key file from all significant sites
def generate_combined_key_file(methylated_dir, output_path):
    combined_keys = set()
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
                    combined_keys.add(key)

    with open(output_path, 'w') as out:
        for key in sorted(combined_keys, key=lambda k: (k.split(':')[0], int(k.split(':')[1]), k.split(':')[2])):
            out.write(key + '\n')

    return combined_keys

# Step 3b: Filter -f files to include only combined keys
def filter_filtered_files_by_keys(filtered_dir, shared_dir, combined_keys):
    os.makedirs(shared_dir, exist_ok=True)

    for filename in os.listdir(filtered_dir):
        if filename.startswith('.'):
            continue

        input_path = os.path.join(filtered_dir, filename)
        output_path = os.path.join(shared_dir, filename)

        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            header = infile.readline()
            outfile.write(header)

            for line in infile:
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    key = f"{fields[0]}:{fields[1]}:{fields[2]}"
                    if key in combined_keys:
                        outfile.write(line)

# Main control
def main():
    parser = argparse.ArgumentParser(description="Filter and merge methylation call files.")
    parser.add_argument("-i", "--input", required=True, help="Input folder path")
    parser.add_argument("-f", "--filtered", required=True, help="Filtered output folder")
    parser.add_argument("-m", "--methylated", required=True, help="Methylated significant output folder")
    parser.add_argument("-s", "--shared", required=True, help="Shared methylation output folder")
    parser.add_argument("-threshold", type=int, required=True, help="Coverage threshold")

    args = parser.parse_args()

    os.makedirs(args.filtered, exist_ok=True)
    os.makedirs(args.methylated, exist_ok=True)
    os.makedirs(args.shared, exist_ok=True)

    print("[1/3] Filtering by total coverage and FDR...")
    for filename in os.listdir(args.input):
        if filename.startswith('.'):
            continue
        input_file = os.path.join(args.input, filename)
        filtered_file = os.path.join(args.filtered, filename)
        filter_file(input_file, filtered_file, args.threshold)

    print("[2/3] Keeping rows with FDR < 0.05...")
    methylation_significant_only(args.filtered, args.methylated)

    print("[3/3] Building combined methylated site list and filtering filtered files...")
    combined_key_file = os.path.join(args.shared, "combined_methylated_sites.txt")
    combined_keys = generate_combined_key_file(args.methylated, combined_key_file)
    filter_filtered_files_by_keys(args.filtered, args.shared, combined_keys)

    print("✅ Done!")

if __name__ == "__main__":
    main()
