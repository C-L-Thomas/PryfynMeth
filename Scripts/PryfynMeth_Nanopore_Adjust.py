#!/usr/bin/env python3

import os
import argparse

def adjust_rows(filepath, output_dir):
    adjusted_lines = []
    with open(filepath, 'r') as file:
        for line in file:
            if not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 6:
                adjusted_lines.append(line)
                continue
            strand = cols[5]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                adjusted_lines.append(line)
                continue

            if strand == '+':
                cols[1] = str(start + 1)
                cols[2] = str(end + 1)
            elif strand == '-':
                cols[1], cols[2] = cols[2], cols[1]

            adjusted_lines.append('\t'.join(cols) + '\n')

    output_path = os.path.join(output_dir, os.path.basename(filepath))
    with open(output_path, 'w') as outfile:
        outfile.writelines(adjusted_lines)

def main():
    parser = argparse.ArgumentParser(description="Adjust nanopore methylation calls.")
    parser.add_argument('-i', '--input_folder', required=True, help='Input folder with dataset files')
    parser.add_argument('-o', '--output_folder', default='adjusted_output', help='Folder to save adjusted files')

    args = parser.parse_args()
    input_folder = args.input_folder
    output_folder = args.output_folder

    os.makedirs(output_folder, exist_ok=True)

    for filename in os.listdir(input_folder):
        if filename.startswith('.'):
            continue  # skip hidden files
        filepath = os.path.join(input_folder, filename)
        if os.path.isfile(filepath):
            adjust_rows(filepath, output_folder)

    print(f"Adjusted files saved in '{output_folder}'.")

if __name__ == '__main__':
    main()
