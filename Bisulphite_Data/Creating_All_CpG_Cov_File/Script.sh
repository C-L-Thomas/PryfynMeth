awk '$6 == "CG" { print $1, $2, $3, 10, 10, $6, $7 }' OFS="\t" file.CpG_report.txt > fake_CpG_report.txt


/Bismark-0.24.2/coverage2cytosine \
  -o master_cpg.cov --merge_CpGs \
  --genome_folder /path_to_genome/  fake_CpG_report.txt 

cut -f1-3 master_cpg.cov.CpG_report.merged_CpG_evidence.cov > destranded_cpgs.cov

