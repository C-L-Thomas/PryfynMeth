import os
import csv
import math
import argparse
from collections import defaultdict

# --- Parse arguments ---
parser = argparse.ArgumentParser(description="Compute methylation means and SEs per condition")
parser.add_argument("-i", "--input", required=True, help="Path to the folder containing input files")
parser.add_argument("-meta", "--metadata", required=True, help="Metadata TSV file with inputFile and condition columns")
parser.add_argument("-o", "--output", default="CpG_condition_means.tsv", help="Output file name")
args = parser.parse_args()

data_folder = args.input
metadata_file = args.metadata
output_file = args.output

# --- Load metadata ---
file_to_condition = {}
with open(metadata_file, 'r') as meta:
    reader = csv.DictReader(meta, delimiter='\t')
    for row in reader:
        file_to_condition[row['inputFile']] = row['condition']

# --- Debug print ---
print("📂 Files in input folder:", os.listdir(data_folder))
print("📋 Files in metadata:", list(file_to_condition.keys()))

# --- Initialize storage ---
methylation_data = defaultdict(lambda: defaultdict(list))
all_positions = set()

# --- Process files ---
for filename in os.listdir(data_folder):
    if filename not in file_to_condition:
        print(f"⚠️ Skipping file not listed in metadata: {filename}")
        continue

    condition = file_to_condition[filename]
    print(f"Processing: {filename} (condition {condition})")

    filepath = os.path.join(data_folder, filename)

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                chr_, pos, strand = row['chr'], int(row['pos']), row['strand']
                C = int(row['C'])
                total = int(row['total'])
            except (ValueError, KeyError):
                continue  # skip malformed rows

            key = (chr_, pos, strand)
            all_positions.add(key)

            if total == 0:
                continue  # exclude from percentage analysis but not from output

            percent = (C / total) * 100
            methylation_data[key][condition].append(percent)

# --- Prepare output ---
all_conditions = sorted(set(file_to_condition.values()), key=int)
output_rows = []

# Header
header = ["chr", "pos", "strand"]
for cond in all_conditions:
    header.extend([f"condition{cond}_mean", f"condition{cond}_SE"])
output_rows.append(header)

# Data rows
for key in sorted(all_positions):
    chr_, pos, strand = key
    row = [chr_, str(pos), strand]
    for cond in all_conditions:
        values = methylation_data[key].get(cond, [])
        if values:
            mean = sum(values) / len(values)
            se = (
                math.sqrt(sum((x - mean) ** 2 for x in values) / len(values)) / math.sqrt(len(values))
                if len(values) > 1 else 0
            )
            row.extend([f"{mean:.2f}", f"{se:.2f}"])
        else:
            row.extend(["NA", "NA"])
    output_rows.append(row)

# --- Write output ---
with open(output_file, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerows(output_rows)

print(f"✅ Done. Output written to {output_file}")
