# First Run dorado (Unless you live call)

dorado basecaller sup,5mCG_5hmCG ./ \
  --reference ref.fa \
  --kit-name SQK-NBD114-24 > calls_sup.bam

# Sort file using Samtools

samtools sort -o calls_sup_sorted.bam calls_sup.bam

# Demultiplex (If using a barcoding kit)
 
dorado demux \
  --output-dir \
  --no-classify \
  calls_sup_sorted.bam

# Prepare cpgs bed file

modkit motif bed reference.fasta CG 0 1> cpgs.bed

# Prepare Nanopore File

python PryfynMeth_Nanopore_Prepare.py -i Input/ -ref cpgs.bed -reduce

# Adjust the output file (as nanopore is 1 off actual genome placement)

python PryfynMeth_Nanopore_Adjust.py -i Methylation_Data/

# Perform Binomial

python PryfynMeth_Binomial.py -meta metadata.txt -platform nano -i Methylation_Data/ -o Methylation_Binomial/

# Filter Results

python PryfynMeth_Filter.py  -i Methylation_Binomial/ -f Filtered_Folder -m Methylated_Folder -s Shared_folder -threshold 10 -revert

# Format for Downstream Differential Analyses

python PryfynMeth_DSS_Prepare.py -i reverted/ -o  DSS/
