import os
import sys
import argparse

def load_reference(ref_file):
    """Loads the reference file into a list of (chrom, start, end, strand) tuples."""
    reference_sites = []
    with open(ref_file, "r") as ref:
        for line in ref:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                reference_sites.append((parts[0], parts[1], parts[2], parts[5]))  # (chrom, start, end, strand)
    return reference_sites

def process_files(input_folder, ref_file):
    """Processes methylation files and ensures all reference sites are included in order."""
    if not os.path.isdir(input_folder):
        print(f"Error: Folder '{input_folder}' not found!")
        sys.exit(1)

    # Load reference CpG sites in order
    reference_sites = load_reference(ref_file)

    # Create output directories in the script's working directory
    script_dir = os.getcwd()
    methylation_folder = os.path.join(script_dir, "Methylation_Data")
    hydroxymethylation_folder = os.path.join(script_dir, "Hydroxymethylation_Data")

    os.makedirs(methylation_folder, exist_ok=True)
    os.makedirs(hydroxymethylation_folder, exist_ok=True)

    # Process each *_meth.bed file
    for filename in os.listdir(input_folder):
        if filename.endswith("_meth.bed"):
            input_file = os.path.join(input_folder, filename)
            prefix = filename.replace("_meth.bed", "")

            print(f"Now processing sample: {prefix}")  # Print sample being processed

            output_h_path = os.path.join(hydroxymethylation_folder, f"{prefix}_CpG_hydroxymeth.bed")
            output_m_path = os.path.join(methylation_folder, f"{prefix}_CpG_meth.bed")

            observed_sites = {}  # Dictionary to store observed sites and their full data
            max_columns = 0  # Track the maximum number of columns found

            with open(input_file, "r") as infile:
                for line in infile:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        chrom, start, end, meth_type = parts[:4]
                        site_key = (chrom, start, end)

                   # Store full line data
                        observed_sites[site_key] = parts

                        # Update max column count if this line has more columns
                        max_columns = max(max_columns, len(parts))

            # Open output files for writing
            with open(output_h_path, "w") as hydrox_file, open(output_m_path, "w") as meth_file:
                for chrom, start, end, strand in reference_sites:
                    site_key = (chrom, start, end)

                    if site_key in observed_sites:
                        # Write observed data
                        full_data = observed_sites[site_key]
                        if full_data[3] == "h":
                            hydrox_file.write("\t".join(full_data) + "\n")
                        elif full_data[3] == "m":
                            meth_file.write("\t".join(full_data) + "\n")
                    else:
                        # Construct missing row with 0s for extra columns
                        missing_columns = ["0"] * (max_columns - 6)  # Ensure correct column count
                        missing_line = [chrom, start, end, ".", ".", strand] + missing_columns
                        hydrox_file.write("\t".join(missing_line) + "\n")
                        meth_file.write("\t".join(missing_line) + "\n")

            print(f"- Created: {output_h_path}")
            print(f"- Created: {output_m_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process methylation bed files and ensure all reference sites are included in order.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input folder containing *_meth.bed files.")
    parser.add_argument("-ref", "--reference", required=True, help="Path to the reference CpG file.")

    args = parser.parse_args()
    process_files(args.input, args.reference)
