import os
import argparse
from collections import defaultdict
from Bio import SeqIO

def find_all_CGs(reference_fasta):
    cg_sites = defaultdict(set)
    print(f"📖 Reading reference genome: {reference_fasta}")
    for record in SeqIO.parse(reference_fasta, "fasta"):
        chrom = record.id.strip()
        seq = str(record.seq).upper()
        for i in range(len(seq) - 1):
            if seq[i:i+2] == "CG":
                cg_sites[chrom].add(i + 1)  # 1-based
    return cg_sites

def load_cov_file(cov_path):
    covered = defaultdict(dict)
    print(f"📄 Loading coverage file: {cov_path}")
    with open(cov_path, "r") as f:
        for line in f:
            if line.strip() == "":
                continue
            parts = line.strip().split("\t")
            if len(parts) != 6:
                print(f"⚠️ Skipping malformed line: {line.strip()}")
                continue
            chrom, start, end, meth_frac, count_meth, count_unmeth = parts
            chrom = chrom.strip()
            try:
                start = int(start)
            except ValueError:
                print(f"⚠️ Skipping non-integer start: {start}")
                continue
            covered[chrom][start] = (meth_frac, count_meth, count_unmeth)
    return covered

def convert_cov_format(cg_sites, covered_sites, output_path):
    with open(output_path, "w") as out:
        written = 0
        for chrom in sorted(cg_sites.keys()):
            if chrom not in covered_sites:
                print(f"[WARN] Chromosome {chrom} not in .cov file — filling with 0s.")
            for pos in sorted(cg_sites[chrom]):
                if pos in covered_sites.get(chrom, {}):
                    _, count_meth, count_unmeth = covered_sites[chrom][pos]
                else:
                    count_meth, count_unmeth = "0", "0"
                out.write(f"{chrom}\t{pos}\t*\t{count_meth}\t{count_unmeth}\tCG\tNA\n")
                written += 1
        print(f"✅ Wrote {written} lines to {output_path}")

def process_all_files(input_dir, output_dir, reference_fasta):
    os.makedirs(output_dir, exist_ok=True)

    print("📖 Indexing CG sites from reference genome...")
    cg_sites = find_all_CGs(reference_fasta)
    total_cg = sum(len(s) for s in cg_sites.values())
    print(f"✅ Found {total_cg} total CG sites across {len(cg_sites)} chromosomes.")
    if total_cg == 0:
        print("❌ No CG sites found! Please check your reference genome.")
        return

    print("📂 Scanning input directory:", input_dir)
    found_cov = False
    for filename in os.listdir(input_dir):
        print(f"🔍 Checking file: {filename}")
        if filename.endswith(".cov"):
            found_cov = True
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename)
            print(f"🔄 Processing {filename}...")
            covered_sites = load_cov_file(input_path)
            print("→ Example .cov chromosomes:", list(covered_sites.keys())[:5])
            convert_cov_format(cg_sites, covered_sites, output_path)

    if not found_cov:
        print("⚠️ No .cov files found in the input directory!")

    print("🏁 Finished processing all files.")

def main():
    parser = argparse.ArgumentParser(description="Convert .cov files to full-CG PryfynMeth format.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing input .cov files")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save converted files")
    parser.add_argument("-ref", "--reference", required=True, help="Reference genome FASTA file")
    args = parser.parse_args()

    process_all_files(args.input_dir, args.output_dir, args.reference)

if __name__ == "__main__":
    main()
