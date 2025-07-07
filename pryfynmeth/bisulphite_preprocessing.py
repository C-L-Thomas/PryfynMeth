import os
import argparse
import pandas as pd

def load_cov_file(cov_path):
    cov_data = set()
    with open(cov_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                cov_data.add(f"{parts[0]}|||{parts[1]}")
    return cov_data

def process_stranded_file(input_path, output_path, output_type):
    modifiers_written = []
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            col1, col2 = parts[0], parts[1]
            strand = parts[2]
            try:
                c_count = int(parts[3])
                t_count = int(parts[4])
            except ValueError:
                c_count = float(parts[3])
                t_count = float(parts[4])
            coverage = c_count + t_count
            col4, col5 = str(c_count), str(coverage)
            col6, col7 = parts[5], parts[6]
            if output_type == 'all':
                outfile.write(f"{col1}\t{col2}\t{strand}\t{col4}\t{col5}\t{col6}\t{col7}\n")
            elif output_type == 'coverage':
                if coverage != 0:
                    outfile.write(f"{col1}\t{col2}\t{strand}\t{col4}\t{col5}\t{col6}\t{col7}\n")
            elif output_type == 'shared':
                if coverage != 0:
                    outfile.write(f"{col1}\t{col2}\t{strand}\t{col4}\t{col5}\t{col6}\t{col7}\n")
                    modifiers_written.append(f"{col1}|||{col2}|||{strand}")
    return modifiers_written

def process_destranded_file(input_path, output_path, output_type, cov_sites=None):
    modifiers_written = []
    input_rows = {}
    with open(input_path, 'r') as infile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            col1 = parts[0]
            col2 = parts[1]
            key = f"{col1}|||{col2}"
            try:
                cov1 = int(parts[4])
                cov2 = int(parts[5])
            except ValueError:
                cov1 = float(parts[4])
                cov2 = float(parts[5])
            total_cov = cov1 + cov2
            row = f"{col1}\t{col2}\t*\t{cov1}\t{total_cov}\t-\t-"
            if output_type == 'all':
                input_rows[key] = row
            elif output_type == 'coverage' and total_cov != 0:
                input_rows[key] = row
            elif output_type == 'shared' and total_cov != 0:
                input_rows[key] = row
                modifiers_written.append(key)

    with open(output_path, 'w') as outfile:
        if output_type == 'all' and cov_sites is not None:
            for key in sorted(cov_sites):
                if key in input_rows:
                    outfile.write(input_rows[key] + '\n')
                else:
                    chrom, pos = key.split('|||')
                    outfile.write(f"{chrom}\t{pos}\t*\t0\t0\t-\t-\n")
        else:
            for row in input_rows.values():
                outfile.write(row + '\n')

    return modifiers_written

def finalize_shared_outputs(output_dir, all_modifiers):
    for filename in os.listdir(output_dir):
        path = os.path.join(output_dir, filename)
        with open(path, 'r') as infile:
            lines = infile.readlines()

        file_mods = set()
        line_dict = {}
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            mod = f"{parts[0]}|||{parts[1]}|||{parts[2]}"
            file_mods.add(mod)
            line_dict[mod] = line.strip()

        final_lines = []
        for mod in sorted(all_modifiers):
            parts = mod.split('|||')
            if len(parts) != 3:
                print(f"⚠️ Skipping malformed modifier: {mod}")
                continue
            col1, col2, col3 = parts
            if mod in line_dict:
                final_lines.append(line_dict[mod])
            else:
                final_lines.append(f"{col1}\t{col2}\t{col3}\t0\t0\t-\t-")

        with open(path, 'w') as outfile:
            for line in final_lines:
                outfile.write(line + '\n')

def main():
    parser = argparse.ArgumentParser(description='PryfynMeth Bisulphite Preprocessing Script')
    parser.add_argument('-i', '--input', required=True, help='Input directory path')
    parser.add_argument('-o', '--output', required=True, help='Output directory path')
    parser.add_argument('-type', required=True, choices=['destranded', 'stranded'], help='Type of input file')
    parser.add_argument('--output_type', required=True, choices=['all', 'shared', 'coverage'], help='Type of output filtering')
    parser.add_argument('cov_file', nargs='?', help='Reference cov file (required if -type destranded and --output_type all)')

    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.output
    file_type = args.type
    output_type = args.output_type
    cov_file = args.cov_file

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if file_type == 'destranded' and output_type == 'all' and not cov_file:
        raise ValueError("You must provide a cov file for -type destranded with --output_type all")

    cov_sites = load_cov_file(cov_file) if file_type == 'destranded' and output_type == 'all' else None

    all_modifiers = set()

    for filename in os.listdir(input_dir):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)

        if file_type == 'stranded':
            modifiers = process_stranded_file(input_path, output_path, output_type)
            if output_type == 'shared':
                all_modifiers.update(modifiers)

        elif file_type == 'destranded':
            modifiers = process_destranded_file(input_path, output_path, output_type, cov_sites)
            if output_type == 'shared':
                all_modifiers.update(modifiers)

    if output_type == 'shared':
        print("Finalizing shared output across all files...")
        finalize_shared_outputs(output_dir, all_modifiers)

if __name__ == '__main__':
    main()
