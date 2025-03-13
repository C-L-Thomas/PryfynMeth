# Background
Welcome to the package PryfynMeth. This aims to be a toolkit to aid with the analysis of insect methylation data. We have included separate subfolders for workthroughs with Whole Genome Bisulphite Sequencing data and Nanopore sequencing data. 

# Preprocessing
For Whole Genome Bisulphite Sequencing files, these scripts function with a particular output of the Bismark Aligner. The specific file type is generated using the Bismark `coverage2cytosine` command, and is a CpG_report file. This is designed for stranded data, but there is no reason it shouldn't work for destranded. An example workflow to generate this file type will be uploaded in future.

For Nanopore sequencing, these scripts function with the outputs of `modkit`. An example usage of modkit will be uploaded. There is a script included that can separate the methylation and hydroxymethylation into different files. This script `PryfynMeth_Nanopore_Separate.py` also adds rows which have 0 coverage in the nanopore bed files. This is important for ensuring the FDR calculation is equal for all samples. However, before running `PryfynMeth_Nanopore_Separate.py`, you must first generate a reference bed for your organism using modkit:

```
modkit motif bed reference.fasta CG 0 1> cg_motifs.bed
```

Once this is generated, you can run `PryfynMeth_Nanopore_Separate.py`:

```
python3 PryfynMeth_Nanopore_Separate.py -i Input/ -ref cg_motifs.bed
```

`-i` - Input folder of modkit output bed files

`-ref` - File generated from modkit motif command outlined above

# Binomial Test 
To determine whether an individual site is methylated in whole genome methylation sequencing, it's common practice to assess if the observed methylation level is significantly higher than what would be expected by chance. This is typically done by comparing the observed methylation proportion to a statistical threshold, often using a binomial test. In whole genome bisulphite sequencing and Nanopore sequencing, this threshold is determined by the percentage methylation found in lambda spiked DNA. To identify methylated sites, you can use the `PryfynMeth_Binomial.py` command. 

```
python3 PryfynMeth_Binomial.py -meta metadata.txt -platform nano -i Path_to_Input_Folder/ -o Path_to_Output_Folder/
```

`-meta` is the input file for your metadata table

`-platform` illu for illumina inputs, nano for nanopore inputs

`-i` the path to a folder containing ONLY the files for binomial processing

`-o` the name of the output folder
