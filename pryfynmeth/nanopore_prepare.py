import os
import sys
import argparse

def load_reference(ref_file):
    reference_sites = []
    with open(ref_file, "r") as ref:
        for line in ref:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                reference_sites.append((parts[0], parts[1], parts[2], parts[5]))  # (chrom, start, end, strand)
    return reference_sites

def extract_modifier(parts):
    if len(parts) >= 6:
        return (parts[0], parts[1], parts[5])  # chrom, start, strand
    return None

def ensure_complete_matrix(folder, suffix, all_modifiers):
    for filename in os.listdir(folder):
        if filename.endswith(f"_{suffix}.bed"):
            path = os.path.join(folder, filename)
            lines_by_mod = {}
            max_cols = 0

            with open(path, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    mod = extract_modifier(parts)
                    if mod:
                        lines_by_mod[mod] = parts
                        max_cols = max(max_cols, len(parts))

            updated_lines = []
            for mod in all_modifiers:
                if mod in lines_by_mod:
                    updated_lines.append("\t".join(lines_by_mod[mod]))
                else:
                    chrom, start, strand = mod
                    end = str(int(start) + 1)
                    missing_row = [chrom, start, end, ".", ".", strand] + ["0"] * (max_cols - 6)
                    updated_lines.append("\t".join(missing_row))

            with open(path, "w") as f:
                f.write("\n".join(updated_lines) + "\n")
            print(f"- Updated with complete matrix: {path}")

def process_files(input_folder, ref_file=None, reduce=False):
    if not os.path.isdir(input_folder):
        print(f"Error: Folder '{input_folder}' not found!")
        sys.exit(1)

    reference_sites = load_reference(ref_file) if ref_file else None

    script_dir = os.getcwd()
    methylation_folder = os.path.join(script_dir, "Methylation_Data")
    hydroxymethylation_folder = os.path.join(script_dir, "Hydroxymethylation_Data")

    os.makedirs(methylation_folder, exist_ok=True)
    os.makedirs(hydroxymethylation_folder, exist_ok=True)

    prefixes = []

    for filename in os.listdir(input_folder):
        if filename.endswith("_meth.bed"):
            input_file = os.path.join(input_folder, filename)
            prefix = filename.replace("_meth.bed", "")
            prefixes.append(prefix)

            print(f"Now processing sample: {prefix}")

            output_h_path = os.path.join(hydroxymethylation_folder, f"{prefix}_CpG_hydroxymeth.bed")
            output_m_path = os.path.join(methylation_folder, f"{prefix}_CpG_meth.bed")

            observed_sites = {}
            max_columns = 0

            with open(input_file, "r") as infile:
                for line in infile:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        chrom, start, end, meth_type = parts[:4]
                        site_key = (chrom, start, end)
                        observed_sites[site_key] = parts
                        max_columns = max(max_columns, len(parts))

            with open(output_h_path, "w") as hydrox_file, open(output_m_path, "w") as meth_file:
                if reference_sites:
                    for chrom, start, end, strand in reference_sites:
                        site_key = (chrom, start, end)
                        if site_key in observed_sites:
                            full_data = observed_sites[site_key]
                            if full_data[3] == "h":
                                hydrox_file.write("\t".join(full_data) + "\n")
                            elif full_data[3] == "m":
                                meth_file.write("\t".join(full_data) + "\n")
                        else:
                            missing_columns = ["0"] * (max_columns - 6)
                            missing_line = [chrom, start, end, ".", ".", strand] + missing_columns
                            hydrox_file.write("\t".join(missing_line) + "\n")
                            meth_file.write("\t".join(missing_line) + "\n")
                else:
                    for parts in observed_sites.values():
                        if parts[3] == "h":
                            hydrox_file.write("\t".join(parts) + "\n")
                        elif parts[3] == "m":
                            meth_file.write("\t".join(parts) + "\n")

            print(f"- Created: {output_h_path}")
            print(f"- Created: {output_m_path}")

    if reduce:
	print("\n-- Generating global modifier lists and enforcing full matrices --")

        for suffix, folder in [("hydroxymeth", hydroxymethylation_folder), ("meth", methylation_folder)]:
            all_modifiers = set()

            # Build global list of unique modifiers
            for filename in os.listdir(folder):
                if filename.endswith(f"_{suffix}.bed"):
                    with open(os.path.join(folder, filename), "r") as f:
                        for line in f:
                            parts = line.strip().split("\t")
                            mod = extract_modifier(parts)
                            if mod:
                                all_modifiers.add(mod)

            all_modifiers = sorted(all_modifiers, key=lambda x: (x[0], int(x[1]), x[2]))  # sort by chrom, start, strand
            ensure_complete_matrix(folder, suffix, all_modifiers)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process methylation bed files with optional reference and matrix reduction.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input folder containing *_meth.bed files.")
    parser.add_argument("-ref", "--reference", required=False, help="(Optional) Path to the reference CpG file.")
    parser.add_argument("-reduce", action="store_true", help="(Optional) Ensure all output files contain identical rows by filling missing sit>

    args = parser.parse_args()
    process_files(args.input, args.reference, args.reduce)
