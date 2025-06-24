import os
import argparse

def process_cov_file(input_path, output_path):
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue  # Skip lines that don't have enough columns

            col1 = parts[0]  # Column 1
            col2 = parts[1]  # Column 2
            star = '*'       # Column 3: '*'
            col4 = parts[4]  # Column 4 becomes original Column 5
            try:
                col5 = str(int(parts[4]) + int(parts[5]))  # Column 5 is Column 5 + Column 6
            except ValueError:
                col5 = str(float(parts[4]) + float(parts[5]))

            extra1 = '-'  # Extra column 1
            extra2 = '-'  # Extra column 2

            outfile.write(f"{col1}\t{col2}\t{star}\t{col4}\t{col5}\t{extra1}\t{extra2}\n")

def process_illu_file(input_path, output_path):
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue  # Skip malformed lines

            col1 = parts[0]  # Chromosome
            col2 = parts[1]  # Position
            star = '*'       # Placeholder
            col4 = parts[3]  # C count
            col5 = str(int(parts[3]) + int(parts[4]))  # Coverage = C + T
            col6 = parts[5]  # Input column 6 (context)
            col7 = parts[6]  # Input column 7 (triplet)

            outfile.write(f"{col1}\t{col2}\t{star}\t{col4}\t{col5}\t{col6}\t{col7}\n")

def main():
    parser = argparse.ArgumentParser(description='PryfynMeth Bisulphite Preprocessing Script')
    parser.add_argument('-i', '--input', required=True, help='Input directory path')
    parser.add_argument('-o', '--output', required=True, help='Output directory path')
    parser.add_argument('-type', required=True, choices=['cov', 'report', 'illu'], help='Type of input file')

    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.output
    file_type = args.type

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)

        if file_type == 'cov':
            process_cov_file(input_path, output_path)
        elif file_type == 'illu':
            process_illu_file(input_path, output_path)
        else:
            print(f"Skipping unsupported type '{file_type}' for file {filename}")

if __name__ == '__main__':
    main()
