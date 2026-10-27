import os
import argparse


# Step 1: Filter based on total coverage and FDR
def filter_file(input_path, output_path, threshold, low_coverage_na=False):

    with open(input_path, "r") as infile, open(output_path, "w") as outfile:

        header = infile.readline()
        outfile.write(header)

        for line in infile:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            try:
                coverage = int(fields[5])
                fdr = float(fields[7])
            except ValueError:
                coverage = 0
                fdr = 1.0

            if coverage < threshold:
                if low_coverage_na:
                    fields[3] = fields[4] = fields[5] = "NA"
                    fields[6] = fields[7] = "NA"
                else:
                    fields[3] = fields[4] = fields[5] = "0"
                    fields[6] = fields[7] = "1.0"

            elif fdr >= 0.05:
                fields[3] = "0"

            outfile.write("\t".join(fields) + "\n")


# Step 2: Keep sites with FDR < 0.05
def methylation_significant_only(filtered_dir, output_dir_m):

    os.makedirs(output_dir_m, exist_ok=True)

    for filename in os.listdir(filtered_dir):

        if filename.startswith("."):
            continue

        input_path = os.path.join(filtered_dir, filename)
        output_path = os.path.join(output_dir_m, filename)

        with open(input_path, "r") as infile, open(output_path, "w") as outfile:

            header = infile.readline()
            outfile.write(header)

            for line in infile:
                fields = line.rstrip("\n").split("\t")

                if len(fields) < 8:
                    continue

                try:
                    fdr = float(fields[7])
                except ValueError:
                    continue

                if fdr < 0.05:
                    outfile.write("\t".join(fields) + "\n")


# Step 3: Build a combined list of sites methylated in at least one sample
def generate_combined_key_file(methylated_dir, output_path):

    combined_keys = set()

    for filename in os.listdir(methylated_dir):

        if filename.startswith("."):
            continue

        path = os.path.join(methylated_dir, filename)

        with open(path, "r") as infile:

            infile.readline()

            for line in infile:
                fields = line.rstrip("\n").split("\t")

                if len(fields) >= 3:
                    key = f"{fields[0]}:{fields[1]}:{fields[2]}"
                    combined_keys.add(key)

    sorted_keys = sorted(
        combined_keys,
        key=lambda k: (
            k.split(":")[0],
            int(k.split(":")[1]),
            k.split(":")[2]
        )
    )

    with open(output_path, "w") as outfile:
        for key in sorted_keys:
            outfile.write(key + "\n")

    return combined_keys


# Step 4: Keep combined sites in every sample
def filter_filtered_files_by_keys(
    filtered_dir,
    shared_dir,
    combined_keys,
    low_coverage_na=False
):

    os.makedirs(shared_dir, exist_ok=True)

    sorted_keys = sorted(
        combined_keys,
        key=lambda k: (
            k.split(":")[0],
            int(k.split(":")[1]),
            k.split(":")[2]
        )
    )

    for filename in os.listdir(filtered_dir):

        if filename.startswith("."):
            continue

        input_path = os.path.join(filtered_dir, filename)
        output_path = os.path.join(shared_dir, filename)

        key_to_line = {}

        with open(input_path, "r") as infile:

            header = infile.readline()

            for line in infile:
                fields = line.rstrip("\n").split("\t")

                if len(fields) >= 3:
                    key = f"{fields[0]}:{fields[1]}:{fields[2]}"
                    key_to_line[key] = line.rstrip("\n")

        with open(output_path, "w") as outfile:

            outfile.write(header)

            for key in sorted_keys:

                if key in key_to_line:
                    outfile.write(key_to_line[key] + "\n")

                else:
                    chrom, pos, strand = key.split(":")

                    if low_coverage_na:
                        outfile.write(
                            f"{chrom}\t{pos}\t{strand}\t"
                            "NA\tNA\tNA\tNA\tNA\n"
                        )
                    else:
                        outfile.write(
                            f"{chrom}\t{pos}\t{strand}\t"
                            "0\t0\t0\t1.0\t1.0\n"
                        )


# Optional step: Restore original counts for adequately covered sites
def final_revert(
    shared_dir,
    input_dir,
    output_dir,
    threshold,
    low_coverage_na=False
):

    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(shared_dir):

        if filename.startswith("."):
            continue

        shared_path = os.path.join(shared_dir, filename)
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)

        # Skip the combined site list and any other files without an input match
        if not os.path.isfile(input_path):
            continue

        input_map = {}

        with open(input_path, "r") as infile:

            infile.readline()

            for line in infile:
                fields = line.rstrip("\n").split("\t")

                if len(fields) >= 7:
                    key = f"{fields[0]}:{fields[1]}:{fields[2]}"
                    input_map[key] = fields

        with open(shared_path, "r") as infile:
            lines = infile.readlines()

        with open(output_path, "w") as outfile:

            outfile.write(lines[0])

            for line in lines[1:]:
                fields = line.rstrip("\n").split("\t")

                if len(fields) >= 7:

                    key = f"{fields[0]}:{fields[1]}:{fields[2]}"

                    try:
                        shared_cov = float(fields[5])
                    except ValueError:
                        shared_cov = 0.0

                    if shared_cov < threshold:

                        missing_value = "NA" if low_coverage_na else "0"

                        fields[3] = missing_value
                        fields[5] = missing_value

                    elif key in input_map:

                        original_fields = input_map[key]

                        fields[3] = original_fields[3]
                        fields[5] = original_fields[5]

                outfile.write("\t".join(fields) + "\n")


def main():

    parser = argparse.ArgumentParser(
        description="Filter and merge methylation call files."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input folder path"
    )

    parser.add_argument(
        "-f", "--filtered",
        required=True,
        help="Filtered output folder"
    )

    parser.add_argument(
        "-m", "--methylated",
        required=True,
        help="Methylated significant output folder"
    )

    parser.add_argument(
        "-s", "--shared",
        required=True,
        help="Shared methylation output folder"
    )

    parser.add_argument(
        "-threshold",
        type=int,
        required=True,
        help="Minimum coverage threshold"
    )

    parser.add_argument(
        "-r", "--reverted",
        help="Reverted output folder; enables the optional revert step"
    )

    parser.add_argument(
        "--low-coverage-na",
        action="store_true",
        help="Use NA for low-coverage and missing sites instead of zero"
    )

    args = parser.parse_args()

    os.makedirs(args.filtered, exist_ok=True)
    os.makedirs(args.methylated, exist_ok=True)
    os.makedirs(args.shared, exist_ok=True)

    print("[1/4] Filtering by total coverage and FDR...")

    for filename in os.listdir(args.input):

        if filename.startswith("."):
            continue

        input_file = os.path.join(args.input, filename)
        filtered_file = os.path.join(args.filtered, filename)

        filter_file(
            input_file,
            filtered_file,
            args.threshold,
            args.low_coverage_na
        )

    print("[2/4] Keeping rows with FDR < 0.05...")

    methylation_significant_only(
        args.filtered,
        args.methylated
    )

    print("[3/4] Building combined methylated site list and shared files...")

    combined_key_file = os.path.join(
        args.shared,
        "combined_methylated_sites.txt"
    )

    combined_keys = generate_combined_key_file(
        args.methylated,
        combined_key_file
    )

    filter_filtered_files_by_keys(
        args.filtered,
        args.shared,
        combined_keys,
        args.low_coverage_na
    )

    if args.reverted:

        print(
            "[4/4] Restoring original counts for sites with "
            f"coverage >= {args.threshold}..."
        )

        final_revert(
            args.shared,
            args.input,
            args.reverted,
            args.threshold,
            args.low_coverage_na
        )

    print("✅ Done!")


if __name__ == "__main__":
    main()
