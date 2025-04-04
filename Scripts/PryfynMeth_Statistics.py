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

            if total > 0:
                covered_rows += 1

            methyl_percent = calculate_methylation_percentage(c, total)
            all_methylation.append(methyl_percent)

            if fdr < 0.05:
                methylated_rows += 1
                methylated_methylation.append(methyl_percent)

    # Averages
    avg_all = sum(all_methylation) / len(all_methylation) if all_methylation else 0.0
    avg_methylated = sum(methylated_methylation) / len(methylated_methylation) if methylated_methylation else 0.0

    # Proportions
    prop_methylated = (methylated_rows / total_rows * 100) if total_rows else 0.0
    prop_covered = (covered_rows / total_rows * 100) if total_rows else 0.0

    return avg_all, avg_methylated, prop_methylated, prop_covered

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
            if filename.endswith("_statistics.txt"):
                continue

            file_path = os.path.join(input_dir, filename)
            if not os.path.isfile(file_path):
                continue

            avg_all, avg_meth, prop_meth, prop_covered = process_file(file_path)

            summary_file.write(f"{filename}\n\n")
            summary_file.write(f"Proportion Methylated (FDR < 0.05): {prop_meth:.2f}%\n")
            summary_file.write(f"Proportion with Sufficient Coverage (total > 0): {prop_covered:.2f}%\n")
            summary_file.write(f"Average Methylation of all sites: {avg_all:.2f}%\n")
            summary_file.write(f"Average Methylation of methylated sites (FDR < 0.05): {avg_meth:.2f}%\n\n")

if __name__ == "__main__":
    main()
