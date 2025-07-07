import os
import argparse

def process_files(input_dir, output_dir):
    """Processes all files in the input directory and writes the output to the output directory."""
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Loop through all files in the input directory
    for file_name in os.listdir(input_dir):
        input_file_path = os.path.join(input_dir, file_name)

        # Skip directories and non-tab-delimited files
        if not os.path.isfile(input_file_path):
            continue

        try:
            # Read the file line by line
            with open(input_file_path, 'r') as infile:
                lines = infile.readlines()

            # Check for a valid header
            if not lines:
                print(f"Skipping file {file_name}: File is empty.")
                continue

            header = lines[0].strip().split('\t')
            required_columns = {"chr", "pos", "total", "C"}

            if not required_columns.issubset(header):
                print(f"Skipping file {file_name}: Missing required columns.")
                continue

            # Get the indices of required columns
            chr_index = header.index("chr")
            pos_index = header.index("pos")
            total_index = header.index("total")
            c_index = header.index("C")

            # Prepare the output file
            output_file_name = f"DSS_{file_name}"
            output_file_path = os.path.join(output_dir, output_file_name)

            with open(output_file_path, 'w') as outfile:
                # Write the new header
                outfile.write("chr\tpos\tN\tX\n")

                # Process and write each line
                for line in lines[1:]:
                    fields = line.strip().split('\t')
                    if len(fields) > max(chr_index, pos_index, total_index, c_index):
                        chr_value = fields[chr_index]
                        pos_value = fields[pos_index]
                        total_value = fields[total_index]
                        c_value = fields[c_index]
                        outfile.write(f"{chr_value}\t{pos_value}\t{total_value}\t{c_value}\n")

            print(f"Processed file: {file_name} -> {output_file_name}")
        except Exception as e:
            print(f"Error processing file {file_name}: {e}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process tab-delimited files for DSS input.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input directory containing the files.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory for the processed files.")

    # Parse arguments
    args = parser.parse_args()

    # Process the files
    process_files(args.input, args.output)
