modkit motif bed reference.fasta CG 0 1> cpgs.bed

python PryfynMeth_Nanopore_Prepare.py -i Input/ -ref cpgs.bed -reduce

python PryfynMeth_Nanopore_Adjust.py -i Methylation_Data/

python PryfynMeth_Binomial.py -meta metadata.txt -platform nano -i Methylation_Data/ -o Methylation_Binomial/

python PryfynMeth_Filter.py  -i Methylation_Binomial/ -f Filtered_Folder -m Methylated_Folder -s Shared_folder -threshold 10 -revert

#python PryfynMeth_DSS_Prepare.py -i reverted/ -o  DSS/
