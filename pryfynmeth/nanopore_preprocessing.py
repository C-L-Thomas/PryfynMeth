import os
import sys
import argparse

# ---------- Helpers ----------

def load_reference(ref_file):
    """
    Load a CpG reference with >=3 columns (chrom, start, end).
    If column 6 exists and is one of '+', '-', '.', use it as strand; else strand='.'.
    Returns a list of (chrom, start, end, strand).
    """
    reference_sites = []
    with open(ref_file, "r") as ref:
        for line in ref:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                chrom, start, end = parts[0], parts[1], parts[2]
                strand = parts[5] if len(parts) >= 6 and parts[5] in ("+", "-", ".") else "."
                reference_sites.append((chrom, start, end, strand))
    print(f"[INFO] Loaded {len(reference_sites):,} reference sites from '{ref_file}'.")
    return reference_sites

def extract_modifier(parts):
    """
    Extract (chrom, start, strand) from a BED-like row with >=6 cols:
      [0]=chrom, [1]=start, [2]=end, [3]=type, [4]=col5, [5]=strand
    """
    if len(parts) >= 6:
        return (parts[0], parts[1], parts[5])
    return None

def fix_col5(parts):
    """Ensure output column 5 (index 4) is never '.', coercing '.' -> '0'."""
    out = parts[:]  # shallow copy
    if len(out) >= 5 and out[4] == ".":
        out[4] = "0"
    return out

# ---------- Core logic ----------

def ensure_complete_matrix(folder, suffix, all_modifiers):
    """
    For each file in 'folder' matching '*_{suffix}.bed', ensure it has exactly the
    rows in 'all_modifiers' (chrom, start, strand), filling missing rows with zeros.
    """
    for filename in os.listdir(folder):
        if not filename.endswith(f"_{suffix}.bed"):
            continue

        path = os.path.join(folder, filename)
        lines_by_mod = {}
        max_cols = 0

        # Read & sanitize existing rows
        with open(path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                parts = fix_col5(parts)
                mod = extract_modifier(parts)
                if mod:
                    lines_by_mod[mod] = parts
                    max_cols = max(max_cols, len(parts))

        updated_lines = []
        if max_cols < 6:
            max_cols = 6  # enforce minimum columns

        # Write in the exact order of all_modifiers, padding where missing
        for mod in all_modifiers:
            if mod in lines_by_mod:
                updated_lines.append("\t".join(lines_by_mod[mod]))
            else:
                chrom, start, strand = mod
                end = str(int(start) + 1)
                # Columns: chrom, start, end, type='.', col5='0', strand, + zero padding
                missing_row = [chrom, start, end, ".", "0", strand] + ["0"] * (max_cols - 6)
                updated_lines.append("\t".join(missing_row))

        with open(path, "w") as f:
            f.write("\n".join(updated_lines) + "\n")
        print(f"- Updated with complete matrix: {path}")

def process_files(input_folder, ref_file=None, reduce=False):
    """
    - If ref_file is provided: use it to drive full-matrix outputs for BOTH meth & hydroxymeth (one line per CpG).
    - If no ref_file: write observed rows only (m -> meth, h -> hydroxymeth).
      (reduce=False means no post-hoc padding; reduce=True can pad both to a master set.)
    """
    if not os.path.isdir(input_folder):
        print(f"Error: Folder '{input_folder}' not found!")
        sys.exit(1)

    reference_sites = load_reference(ref_file) if ref_file else None
    if ref_file and (not reference_sites):
        raise ValueError(
            f"Reference provided ('{ref_file}') but 0 sites were loaded. "
            "Ensure it's tab-delimited with at least 3 columns (chrom, start, end)."
        )

    script_dir = os.getcwd()
    methylation_folder = os.path.join(script_dir, "Methylation_Data")
    hydroxymethylation_folder = os.path.join(script_dir, "Hydroxymethylation_Data")

    os.makedirs(methylation_folder, exist_ok=True)
    os.makedirs(hydroxymethylation_folder, exist_ok=True)

    for filename in os.listdir(input_folder):
        if not filename.endswith("_meth.bed"):
            continue

        input_file = os.path.join(input_folder, filename)
        prefix = filename[:-9]  # strip "_meth.bed"

        print(f"Now processing sample: {prefix}")

        output_h_path = os.path.join(hydroxymethylation_folder, f"{prefix}_CpG_hydroxymeth.bed")
        output_m_path = os.path.join(methylation_folder, f"{prefix}_CpG_meth.bed")

        # Keep *separate* maps so h and m at same site don't overwrite each other
        observed_h = {}
        observed_m = {}
        max_columns = 0

        # Read input once, sanitize col5
        with open(input_file, "r") as infile:
            for line in infile:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    parts = fix_col5(parts)
                    chrom, start, end, meth_type = parts[:4]
                    site_key = (chrom, start, end)
                    max_columns = max(max_columns, len(parts))
                    if meth_type == "h":
                        observed_h[site_key] = parts
                    elif meth_type == "m":
                        observed_m[site_key] = parts
                    else:
                        # ignore unknown types
                        pass

        with open(output_h_path, "w") as hydrox_file, open(output_m_path, "w") as meth_file:
            if reference_sites is not None:
                # Use reference to produce full matrices for BOTH outputs (one line per CpG each)
                base_pad_len = max(6, max_columns) - 6
                for chrom, start, end, strand in reference_sites:
                    site_key = (chrom, start, end)
                    pad = ["0"] * base_pad_len
                    zero_line = [chrom, start, end, ".", "0", strand] + pad

                    # Hydroxy line
                    if site_key in observed_h:
                        hydrox_file.write("\t".join(observed_h[site_key]) + "\n")
                    else:
                        hydrox_file.write("\t".join(zero_line) + "\n")

                    # Methyl line
                    if site_key in observed_m:
                        meth_file.write("\t".join(observed_m[site_key]) + "\n")
                    else:
                        meth_file.write("\t".join(zero_line) + "\n")
            else:
                # No reference: just split observed rows
                for parts in observed_h.values():
                    hydrox_file.write("\t".join(parts) + "\n")
                for parts in observed_m.values():
                    meth_file.write("\t".join(parts) + "\n")

        print(f"- Created: {output_h_path}")
        print(f"- Created: {output_m_path}")

    # Optional: post-hoc normalization to a single master set across both folders
    if reduce:
        print("\n-- Generating master modifier list and enforcing full matrices in BOTH folders --")

        # Build master modifiers from the reference if provided; else from methylation outputs
        if reference_sites:
            master_modifiers = [(c, s, st) for (c, s, e, st) in reference_sites]
        else:
            # derive from methylation outputs
            master_modifiers = set()
            for fname in os.listdir(methylation_folder):
                if not fname.endswith("_CpG_meth.bed"):
                    continue
                with open(os.path.join(methylation_folder, fname), "r") as f:
                    for line in f:
                        parts = line.strip().split("\t")
                        parts = fix_col5(parts)
                        mod = extract_modifier(parts)
                        if mod:
                            master_modifiers.add(mod)
            master_modifiers = sorted(master_modifiers, key=lambda x: (x[0], int(x[1]), x[2]))

        ensure_complete_matrix(methylation_folder, "meth", master_modifiers)
        ensure_complete_matrix(hydroxymethylation_folder, "hydroxymeth", master_modifiers)

# ---------- CLI ----------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process methylation bed files with optional reference-driven full matrices."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Path to the input folder containing *_meth.bed files.")
    parser.add_argument("-ref", "--reference", required=False,
                        help="(Optional) Path to the reference CpG file (tab-delimited; >=3 cols).")
    parser.add_argument("-reduce", action="store_true",
                        help="(Optional) After writing outputs, enforce identical rows across ALL files in BOTH folders.")

    args = parser.parse_args()
    process_files(args.input, args.reference, args.reduce)
