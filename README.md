# Background
Welcome to the package PryfynMeth. This aims to be a toolkit to aid with the analysis of insect methylation data. We have included separate subfolders for workthroughs with Whole Genome Bisulphite Sequencing data and Nanopore sequencing data, and example scripts for the whole pipeline for [Nanopore](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Nanopore_WorkFlow/Example_Workflow.sh).

# Preprocessing
PryfynMeth can accept three input file types; Nanopore [MethylBed](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Nanopore_WorkFlow/example_methyl.bed), [Bismark bismark_methylation_extractor reports (bisulphite stranded)](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Bisulphite_Data/example_stranded.report.txt), and [Bismark coverage2cytosine cov files (bisulphite destranded)](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Bisulphite_Data/example_destranded.cov). 

## Preprocessing: Whole Genome Bisulphite Sequencing

For Whole Genome Bisulphite Sequencing files, these scripts function with a particular output of the Bismark Aligner. There are two file type options; stranded and destranded. For stranded files, the file type is generated using Bismarks bismark_methylation_extractor command with the --report option. This will output a file that ends in *report.txt and is in [this](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Bisulphite_Data/example_stranded.report.txt) format. This is the file input format for the binomial test `PryfynMeth_Binomial.py`.

If you wish to examine destranded data, the input file is different, and there will be an exra processing step to get the file in a similar format to stranded data. Destranded file types are *.cov output files from Bismark's coverage2cytosine command with the --merge_CpGs option used. This will generate files in [this](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Bisulphite_Data/example_destranded.cov) format. Once you have generated cov files, run the following command:

```
python PryfynMeth_Bisulphite_Preprocessing.py -i Input -o all_Input -type stranded --output_type all #template_cov_file.cov

```
`-i` - Input folder cov files in [this](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Bisulphite_Data/example_destranded.cov) format or report files in [this](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Bisulphite_Data/example_stranded.report.txt) format.

`-o` - Output folder for reformatted files. This will be your input folder for `PryfynMeth_Binomial.py`

`-type` - Stranded or Destranded

`-output_type` all, coverage or shared. Explained [here](). Note, if using destranded data, you will need an input "template" cov file in the position indicated by the hash. How to generate this will be explained [here]().

## Preprocessing: Nanopore Sequencing

For Nanopore sequencing, these scripts function with the outputs of `modkit`. An example usage of modkit will be uploaded. There is a script included that can separate the methylation and hydroxymethylation into different files. This script `PryfynMeth_Nanopore_Separate.py` also adds rows which have 0 coverage in the nanopore bed files. This is important for ensuring the FDR calculation is equal for all samples. However, before running `PryfynMeth_Nanopore_Separate.py`, you must first generate a reference bed for your organism using modkit:

```
modkit motif bed reference.fasta CG 0 1> cg_motifs.bed
```

As a default, Modkit outputs will give you methylation and hydroxymethylation for sites with coverage. This leads to the problem that you are essentially doubling the number of input rows, which will have an impact on the binomial FDR. The `PryfynMeth_Nanopore_Preprepare.py` command separates these into separate files, thus reducing the input files for FDR. There are three ways you can then decide to deal with this. The first is to proceed with with only the sites with coverage. This requires you only split methylation from hydroxymethylation. To run this:

```
python3 PryfynMeth_Nanopore_Prepare.py -i Input/
````

You may decide that this biases samples with a lower number of sites in the FDR. To resolve this, you may take one of two approaches. The first is to make sure that the files that will go through the binomial possess every CpG in the genome:

```
python3 PryfynMeth_Nanopore_Prepare.py -i Input/ -ref cpgs.bed 
````

Whilst this makes sure each sample is treated the same in the FDR, it does increase the number of sites included in the FDR calculation which may result in fewer sites being classified as methylated. To combat this whilst treating each sample equally, you may wish to remove genomic sites that don't posses coverage in any sample using the -reduce option:

```
python3 PryfynMeth_Nanopore_Prepare.py -i Input/ -ref cpgs.bed -reduce
````

`-i` - Input folder of modkit output bed files, which are in [this](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Nanopore_WorkFlow/example_methyl.bed) format

`-ref` - File generated from modkit motif command outlined above

`-reduce` - Keeps same number of rows between all samples, but removes sites that never have coverage

A final step for the Nanopore preprocessing is to adjust the CG positions. For an unknown reason (I have contacted ONT), modkit files come out 1 position out of sync for + strand bases, and the start and end position is the wrong way around for - strand bases. The PryfynMeth_Nanopore_Adjust.py command fixes this before the binomial test:

```
python PryfynMeth_Nanopore_Adjust.py -i Input/
````

`-i` - Input folder (output from Nanopore_Prepare)

This should generate a folder named adjusted_output.

# Binomial Test 
To determine whether an individual site is methylated in whole genome methylation sequencing, it's common practice to assess if the observed methylation level is significantly higher than what would be expected by chance. This is typically done by comparing the observed methylation proportion to a statistical threshold, often using a binomial test. In whole genome bisulphite sequencing and Nanopore sequencing, this threshold is determined by the percentage methylation found in lambda spiked DNA. To identify methylated sites, you can use the `PryfynMeth_Binomial.py` command. 

```
python3 PryfynMeth_Binomial.py -meta metadata.txt -platform nano -i Path_to_Input_Folder/ -o Path_to_Output_Folder/
```

`-meta` is the input file for your metadata table. An example of the input metadata table can be found [here](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Nanopore_WorkFlow/Example_Nanopore_Metadata.txt)

`-platform` illu for illumina inputs, nano for nanopore inputs

`-i` the path to a folder containing ONLY the files for binomial processing

`-o` the name of the output folder

# Filtering
Once you have generated binomial results, you may wish to filter samples with low read count. The PryfynMeth_Filter.py command does just that. By setting a threshold it will output three folders. The first will be a full list of your sites in your binomial test output, but with the values and statistics adjusted considering your filtering. The second will be a folde of each sample's methylated sites (sites with FDR less than 0.05). The final, will scan each of your methylated site files, and will add genomic locations that have been excluded from other samples, making sure each sample has the same number of input site.

```
python PryfynMeth_Filter.py  -i Binomial_Results_Folder -f output_adjusted_sites -m output_pass_binomial -s output_shared_sites -threshold 10 -revert

```

`-i` An input folder which should be the output (-o) from the binomial step

`-f` An output folder that will give a full list of sites, but with adjusted values. For example, if threshold is set at 10 and the site in question passes FDR but has fewer reads than 10, the reads, methylation count and non-methylation could will be set at 0, and the FDR will become 1 (failing). 

`-m` An output folder that subsets each input file to only give methylated sites (ie those with an FDR less than 0.05).

`-threshold` Anything lower than the number set here will convert the methylation count, non-methylation count and coverage numbers to 0, making the FDR fail.

`-s` An output folder that subsets the input files, so that only rows that are methylated in at least one sample are kept. 

`-revert` This is an optional tag. If included, sites which failed binomial but have over the threshold coverage will be reverted to their original input values. The logic for this is that by changing these values to 0, you may artificially be increasing the difference between samples and consequentially increase the number of differentially methylated sites. This example has threshold set at 10: 

| Input  | Default | -revert |
| ------------- | ------------- | -|
| 1 / 9 (fail FDR) | 0/0  | 0/0 |
| 1/20 (fail FDR)  | 0/0  | 1/20 |
| 2/20 (pass FDR)  | 2/20  | 2/20 |

I find it easiest to run all options (filtered, methylated, subset and revert), which allows you to inspect how many sites proceed for each option, before making a decision for the next step.

# Statistics
To obtain the percentage of methylated sites, level of methylation and the coverage, you can use PryfynMeth_Statistics.py. The input folder (-i) can be your binomial test results, your filtered results, methylated sites results or shared CpG results (although the latter 2 won't give you measures as to how much of the unmethylated genome you're missing).

```
python3 PryfynMeth_Statistics.py -i input_folder

```

# Prepare For Differential Methylation Analysis

For use with [DSS](https://www.bioconductor.org/packages/devel/bioc/vignettes/DSS/inst/doc/DSS.html), you can use the PryfynMeth_DSS_Prepare.py command on either your binomial, filtered, methylated or shared folders:

```
python3 PryfynMeth_DSS_Prepare.py -i output_shared_sites -o DSS_folder

```

`-i` An input folder with the desired files to be used for DSS (eg. -f, -m or -s output from PryfynMeth_Filter.py)

`-o` An output folder for files in the DSS input format

# Additional Tools

PryfynMeth_Condition_CpG_Mean.py takes an input [metadata](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Condition_Summaries/Condition_Metadata.txt) file with your sample conditions, and gives a mean (and SE) for each sites methylation. The output is a file with each CpG with each conditions mean and SE. It is advised the -f folder from PryfynMeth_Filter.py is used as the input.

```
python3 PryfynMeth_Condition_CpG_Mean.py -i input_folder -meta condition_metadata.txt -o filename.tsv -g gene_info.txt
```

PryfynMeth_Sample_Combine.py can take your preprocessed data with an input [metadata](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Condition_Summaries/combine_metadata.txt) to let you determine methylated sites across all samples. This should be proceeded with a binmoial test and filtering. 

```
PryfynMeth_Sample_Combine.py -i Preprocessed_Data/ -m combine_metadata.txt -o Output_Folder/
```
