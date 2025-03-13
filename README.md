# Background
Welcome to the package PryfynMeth. This aims to be a toolkit to aid with the analysis of insect methylation data. We have included separate subfolders for workthroughs with Whole Genome Bisulphite Sequencing data and Nanopore sequencing data. 

# Binomial Test 
To determine whether an individual site is methylated in whole genome methylation sequencing, it's common practice to assess if the observed methylation level is significantly higher than what would be expected by chance. This is typically done by comparing the observed methylation proportion to a statistical threshold, often using a binomial test. In whole genome bisulphite sequencing and Nanopore sequencing, this threshold is determined by the percentage methylation found in lambda spiked DNA. To identify methylated sites, you can use the `PryfynMeth_Binomial.py` command. 

```
python3 PryfynMeth_Binomial.py -meta metadata.txt -platform nano -i Path_to_Input_Folder/ -o Path_to_Output_Folder/
```

`-meta` is the input file for your metadata table

`-platform` illu for illumina inputs, nano for nanopore inputs

`-i` the path to a folder containing ONLY the files for binomial processing

`-o` the name of the output folder
