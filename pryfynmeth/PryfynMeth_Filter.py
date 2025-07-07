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
                fdr = float(fields[7])
            except ValueError:
                coverage = 0
                fdr = 1.0

            if coverage < threshold:
                fields[3] = fields[4] = fields[5] = '0'
                fields[6] = fields[7] = '1.0'
            elif fdr >= 0.05:
                fields[3] = '0'

            outfile.write('\t'.join(fields) + '\n')

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

def filter_filtered_files_by_keys(filtered_dir, shared_dir, combined_keys):
    os.makedirs(shared_dir, exist_ok=True)

    for filename in os.listdir(filtered_dir):
        if filename.startswith('.'):
            continue

        input_path = os.path.join(filtered_dir, filename)
        output_path = os.path.join(shared_dir, filename)

        key_to_line = {}
        with open(input_path, 'r') as infile:
            header = infile.readline()
            for line in infile:
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    key = f"{fields[0]}:{fields[1]}:{fields[2]}"
                    key_to_line[key] = line.strip()

        with open(output_path, 'w') as outfile:
            outfile.write(header)
            for key in sorted(combined_keys, key=lambda k: (k.split(':')[0], int(k.split(':')[1]), k.split(':')[2])):
                if key in key_to_line:
                    outfile.write(key_to_line[key] + '\n')
                else:
                    chrom, pos, strand = key.split(':')
                    outfile.write(f"{chrom}\t{pos}\t{strand}\t0\t0\t0\t1.0\t1.0\n")

def final_revert(shared_dir, input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(shared_dir):
        if filename.startswith('.'):
            continue

        shared_path = os.path.join(shared_dir, filename)
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)

        if not os.path.exists(input_path):
            continue

        input_map = {}
        with open(input_path, 'r') as f:
            _ = f.readline()
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 7:
                    key = f"{fields[0]}:{fields[1]}:{fields[2]}"
                    input_map[key] = fields

        with open(shared_path, 'r') as f:
            lines = f.readlines()

        with open(output_path, 'w') as f:
            f.write(lines[0])
            for line in lines[1:]:
                fields = line.strip().split('\t')
                if len(fields) >= 7:
                    key = f"{fields[0]}:{fields[1]}:{fields[2]}"
                    try:
                        shared_cov = float(fields[5])
                    except ValueError:
                        shared_cov = 0.0

                    if shared_cov < 10:
                        fields[3] = '0'
                        fields[5] = '0'
                    elif key in input_map:
                        i_fields = input_map[key]
                        if fields[3] != i_fields[3]:
                            fields[3] = i_fields[3]
                        if fields[5] != i_fields[5]:
                            fields[5] = i_fields[5]

                f.write('\t'.join(fields) + '\n')

def main():
    parser = argparse.ArgumentParser(description="Filter and merge methylation call files.")
    parser.add_argument("-i", "--input", required=True, help="Input folder path")
    parser.add_argument("-f", "--filtered", required=True, help="Filtered output folder")
    parser.add_argument("-m", "--methylated", required=True, help="Methylated significant output folder")
    parser.add_argument("-s", "--shared", required=True, help="Shared methylation output folder")
    parser.add_argument("-threshold", type=int, required=True, help="Coverage threshold")
    parser.add_argument("-revert", action="store_true", help="Apply final revert step using original input")

    args = parser.parse_args()

    os.makedirs(args.filtered, exist_ok=True)
    os.makedirs(args.methylated, exist_ok=True)
    os.makedirs(args.shared, exist_ok=True)

    print("[1/4] Filtering by total coverage and FDR...")
    for filename in os.listdir(args.input):
        if filename.startswith('.'):
            continue
        input_file = os.path.join(args.input, filename)
        filtered_file = os.path.join(args.filtered, filename)
        filter_file(input_file, filtered_file, args.threshold)

    print("[2/4] Keeping rows with FDR < 0.05...")
    methylation_significant_only(args.filtered, args.methylated)

    print("[3/4] Building combined methylated site list and filtering filtered files...")
    combined_key_file = os.path.join(args.shared, "combined_methylated_sites.txt")
    combined_keys = generate_combined_key_file(args.methylated, combined_key_file)
    filter_filtered_files_by_keys(args.filtered, args.shared, combined_keys)

    if args.revert:
        print("[4/4] Reverting column 4 and 6 from input if values differ (with coverage >= 10 condition)...")
        reverted_dir = os.path.join(os.path.dirname(args.shared), "reverted")
        final_revert(args.shared, args.input, reverted_dir)

    print("\u2705 Done!")

if __name__ == "__main__":
    main()
