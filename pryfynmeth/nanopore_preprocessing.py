import os
import sys
import argparse
from collections import defaultdict

def load_reference(ref_file):
    """
    Read a reference BED-like file and return a list of tuples:
    (chrom, start, end, strand)
    Requires >= 6 columns so strand is parts[5].
    """
    reference_sites = []
    with open(ref_file, "r") as ref:
        for line in ref:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6:
                reference_sites.append((parts[0], parts[1], parts[2], parts[5]))
    return reference_sites

def extract_modifier(parts):
    """
    For matrix harmonization: pull (chrom, start, strand).
    Only valid if the row has at least 6 columns.
    """
    if len(parts) >= 6:
        return (parts[0], parts[1], parts[5])  # chrom, start, strand
    return None

def ensure_complete_matrix(folder, suffix, all_modifiers):
    """
    For every file in `folder` ending with _{suffix}.bed,
    ensure it has exactly one row for each modifier in `all_modifiers`.
    Missing rows are filled with placeholders having correct width.
    """
    for filename in os.listdir(folder):
        if filename.endswith(f"_{suffix}.bed"):
            path = os.path.join(folder, filename)
            lines_by_mod = {}
            max_cols = 0

            with open(path, "r") as f:
                for line in f:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    mod = extract_modifier(parts)
                    if mod:
                        # keep the first occurrence; if duplicates exist, last one wins
                        lines_by_mod[mod] = parts
                        max_cols = max(max_cols, len(parts))

            updated_lines = []
            # Non-negative guard for column padding
            pad = max(0, max_cols - 6)

            for mod in all_modifiers:
                if mod in lines_by_mod:
                    updated_lines.append("\t".join(lines_by_mod[mod]))
                else:
                    chrom, start, strand = mod
                    end = str(int(start) + 1)
                    missing_row = [chrom, start, end, ".", ".", strand] + ["0"] * pad
                    updated_lines.append("\t".join(missing_row))

            with open(path, "w") as f:
                f.write("\n".join(updated_lines) + ("\n" if updated_lines else ""))
            print(f"- Updated with complete matrix: {path}")

def process_files(input_folder, ref_file=None, reduce=False):
    if not os.path.isdir(input_folder):
        print(f"Error: Folder '{input_folder}' not found!")
        sys.exit(1)

    reference_sites = load_reference(ref_file) if ref_file else None

    # Use current working directory for output subfolders
    script_dir = os.getcwd()
    methylation_folder = os.path.join(script_dir, "Methylation_Data")
    hydroxymethylation_folder = os.path.join(script_dir, "Hydroxymethylation_Data")

    os.makedirs(methylation_folder, exist_ok=True)
    os.makedirs(hydroxymethylation_folder, exist_ok=True)

    for filename in sorted(os.listdir(input_folder)):
        if filename.endswith("_meth.bed"):
            input_file = os.path.join(input_folder, filename)
            prefix = filename[:-9]  # drop "_meth.bed"

            print(f"Now processing sample: {prefix}")

            output_h_path = os.path.join(hydroxymethylation_folder, f"{prefix}_CpG_hydroxymeth.bed")
            output_m_path = os.path.join(methylation_folder, f"{prefix}_CpG_meth.bed")

            # Keep both modifiers per (chrom,start,end)
            observed_sites = defaultdict(dict)  # site_key -> {'h': parts, 'm': parts}
            max_columns = 0

            with open(input_file, "r") as infile:
                for line in infile:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 4:
                        chrom, start, end, meth_type = parts[:4]
                        site_key = (chrom, start, end)
                        if meth_type in ("h", "m"):
                            observed_sites[site_key][meth_type] = parts
                        max_columns = max(max_columns, len(parts))

            pad = max(0, max_columns - 6)

            with open(output_h_path, "w") as hydrox_file, open(output_m_path, "w") as meth_file:
                if reference_sites:
                    # Emit exactly one row per reference site, per output
                    for chrom, start, end, strand in reference_sites:
                        site_key = (chrom, start, end)

                        # Hydroxymethylation
                        if site_key in observed_sites and 'h' in observed_sites[site_key]:
                            hydrox_file.write("\t".join(observed_sites[site_key]['h']) + "\n")
                        else:
                            missing_line_h = [chrom, start, end, ".", ".", strand] + ["0"] * pad
                            hydrox_file.write("\t".join(missing_line_h) + "\n")

                        # Methylation
                        if site_key in observed_sites and 'm' in observed_sites[site_key]:
                            meth_file.write("\t".join(observed_sites[site_key]['m']) + "\n")
                        else:
                            missing_line_m = [chrom, start, end, ".", ".", strand] + ["0"] * pad
                            meth_file.write("\t".join(missing_line_m) + "\n")
                else:
                    # No reference: write whatever exists, split by modifier
                    # (Order will be arbitrary unless you want to sort)
                    for per_site in observed_sites.values():
                        if 'h' in per_site:
                            hydrox_file.write("\t".join(per_site['h']) + "\n")
                        if 'm' in per_site:
                            meth_file.write("\t".join(per_site['m']) + "\n")

            print(f"- Created: {output_h_path}")
            print(f"- Created: {output_m_path}")

    if reduce:
        print("\n-- Generating global modifier lists and enforcing full matrices --")
        for suffix, folder in [("hydroxymeth", hydroxymethylation_folder), ("meth", methylation_folder)]:
            all_modifiers = set()

            # Build global list of unique modifiers across all sample files
            for fname in os.listdir(folder):
                if fname.endswith(f"_{suffix}.bed"):
                    with open(os.path.join(folder, fname), "r") as f:
                        for line in f:
                            if not line.strip() or line.startswith("#"):
                                continue
                            parts = line.rstrip("\n").split("\t")
                            mod = extract_modifier(parts)
                            if mod:
                                all_modifiers.add(mod)

            # Sort by chrom, start (numeric), then strand
            all_modifiers = sorted(all_modifiers, key=lambda x: (x[0], int(x[1]), x[2]))
            ensure_complete_matrix(folder, suffix, all_modifiers)

def main():
    parser = argparse.ArgumentParser(
        description="Split *_meth.bed into CpG_hydroxymeth and CpG_meth with optional reference and matrix reduction."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the input folder containing *_meth.bed files.")
    parser.add_argument("-ref", "--reference", required=False, help="(Optional) Path to the reference CpG file (>=6 columns; strand in col 6).")
    parser.add_argument("-reduce", action="store_true", help="(Optional) Ensure all output files contain identical rows by filling missing sites.")
    args = parser.parse_args()

    process_files(args.input, args.reference, args.reduce)

if __name__ == "__main__":
    main()
