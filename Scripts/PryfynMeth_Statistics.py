#!/usr/bin/env python3

import os
import argparse

def calculate_methylation_percentage(c, total):
    if total == 0:
        return 0.0
    return (c / total) * 100

def process_file(file_path):
    all_methylation = []
    methylated_methylation = []

    total_rows = 0
    methylated_rows = 0
    covered_rows = 0
    cov10 = cov20 = cov30 = cov50 = cov100 = 0

    total_sum = 0
    valid_total_rows = 0

    with open(file_path, 'r') as infile:
        header = infile.readline()
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            try:
                c = int(parts[3])
                total = int(parts[5])
                fdr = float(parts[7])
            except ValueError:
                continue

            total_rows += 1
            total_sum += total
            valid_total_rows += 1

            if total > 0:
                covered_rows += 1
            if total >= 10:
                cov10 += 1
            if total >= 20:
                cov20 += 1
            if total >= 30:
                cov30 += 1
            if total >= 50:
                cov50 += 1
            if total >= 100:
                cov100 += 1

            methyl_percent = calculate_methylation_percentage(c, total)
            all_methylation.append(methyl_percent)

            if fdr < 0.05:
                methylated_rows += 1
                methylated_methylation.append(methyl_percent)

    def percent(n): return (n / total_rows * 100) if total_rows else 0.0

    return {
        "avg_all": sum(all_methylation) / len(all_methylation) if all_methylation else 0.0,
        "avg_meth": sum(methylated_methylation) / len(methylated_methylation) if methylated_methylation else 0.0,
        "prop_meth": percent(methylated_rows),
        "prop_cov": percent(covered_rows),
        "cov10": percent(cov10),
        "cov20": percent(cov20),
        "cov30": percent(cov30),
        "cov50": percent(cov50),
        "cov100": percent(cov100),
        "coverage": (total_sum / valid_total_rows) if valid_total_rows else 0.0
    }

def main():
    parser = argparse.ArgumentParser(description="Calculate methylation statistics from input files.")
    parser.add_argument('-i', '--input', required=True, help="Input folder containing data files.")
    args = parser.parse_args()

    input_dir = args.input
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a directory.")
        return

    summary_path = os.path.join(input_dir, "statistics.txt")
    with open(summary_path, 'w') as summary_file:
        for filename in sorted(os.listdir(input_dir)):
            if filename.endswith("_statistics.txt") or filename == "statistics.txt":
                continue

            file_path = os.path.join(input_dir, filename)
            if not os.path.isfile(file_path):
                continue

            stats = process_file(file_path)

            summary_file.write(f"{filename}\n\n")
            summary_file.write(f"Proportion Methylated (FDR < 0.05): {stats['prop_meth']:.2f}%\n")
            summary_file.write(f"Proportion with Sufficient Coverage (total > 0): {stats['prop_cov']:.2f}%\n")
            summary_file.write(f"Proportion with ≥10x coverage: {stats['cov10']:.2f}%\n")
            summary_file.write(f"Proportion with ≥20x coverage: {stats['cov20']:.2f}%\n")
            summary_file.write(f"Proportion with ≥30x coverage: {stats['cov30']:.2f}%\n")
            summary_file.write(f"Proportion with ≥50x coverage: {stats['cov50']:.2f}%\n")
            summary_file.write(f"Proportion with ≥100x coverage: {stats['cov100']:.2f}%\n")
            summary_file.write(f"Average Methylation of all sites: {stats['avg_all']:.2f}%\n")
            summary_file.write(f"Average Methylation of methylated sites (FDR < 0.05): {stats['avg_meth']:.2f}%\n")
            summary_file.write(f"Average Coverage: {stats['coverage']:.2f}%\n\n")

if __name__ == "__main__":
    main()
