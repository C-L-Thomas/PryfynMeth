#!/usr/bin/env python3
import argparse
import os
import glob
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="Build methylated & coverage count matrices from *_binomial_results.txt files"
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Directory containing *_binomial_results.txt files"
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="Output directory to write methylated_counts_matrix.tsv and coverage_counts_matrix.tsv"
    )
    args = parser.parse_args()

    pattern = os.path.join(args.input_dir, "*_binomial_results.txt")
    files = sorted(glob.glob(pattern))

    if not files:
        raise SystemExit(f"No files matching *_binomial_results.txt in {args.input_dir}")

    os.makedirs(args.output_dir, exist_ok=True)

    methylated_series_list = []
    coverage_series_list = []

    for path in files:
        fname = os.path.basename(path)
        sample = fname.replace("_binomial_results.txt", "")

        df = pd.read_csv(path, sep="\t")

        # Basic required columns
        for col in ("chr", "pos"):
            if col not in df.columns:
                raise SystemExit(f"{path}: missing required column '{col}'. Found: {list(df.columns)}")

        # Strand may be present in some formats; fall back to '.' if absent
        strand_col = "strand" if "strand" in df.columns else None

        if strand_col:
            cpg_id = df["chr"].astype(str) + "_" + df["pos"].astype(str) + "_" + df["strand"].astype(str)
        else:
            cpg_id = df["chr"].astype(str) + "_" + df["pos"].astype(str)

        # --- Methylated counts ---
        # Prefer X (DSS-style), else C, else 5th column
        if "X" in df.columns:
            meth = df["X"]
        elif "C" in df.columns:
            meth = df["C"]
        else:
            if df.shape[1] < 5:
                raise SystemExit(f"{path}: need at least 5 columns to take methylated counts by position.")
            meth = df.iloc[:, 4]  # 5th column (0-based index 4)

        s_meth = pd.to_numeric(meth, errors="coerce").fillna(0).astype(int)
        s_meth.index = cpg_id
        s_meth.name = sample
        methylated_series_list.append(s_meth)

        # --- Coverage counts ---
        # Prefer N (DSS-style), else take 6th column as requested
        if "N" in df.columns:
            cov = df["N"]
        else:
            if df.shape[1] < 6:
                raise SystemExit(f"{path}: need at least 6 columns to take coverage by position (6th column).")
            cov = df.iloc[:, 5]  # 6th column (0-based index 5)

        s_cov = pd.to_numeric(cov, errors="coerce").fillna(0).astype(int)
        s_cov.index = cpg_id
        s_cov.name = sample
        coverage_series_list.append(s_cov)

    # Combine into matrices (outer join on CpG IDs)
    meth_mat = pd.concat(methylated_series_list, axis=1).fillna(0).astype(int)
    cov_mat  = pd.concat(coverage_series_list, axis=1).fillna(0).astype(int)

    meth_mat.index.name = "CpG"
    cov_mat.index.name  = "CpG"

    # Write outputs
    meth_out = os.path.join(args.output_dir, "methylated_counts_matrix.tsv")
    cov_out  = os.path.join(args.output_dir, "coverage_counts_matrix.tsv")

    meth_mat.to_csv(meth_out, sep="\t")
    cov_mat.to_csv(cov_out, sep="\t")

if __name__ == "__main__":
    main()
