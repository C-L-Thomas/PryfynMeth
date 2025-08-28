# Background
Welcome to the package PryfynMeth. PryfynMeth is a combination of tools designed to help the analysis of insect methylation data. It takes aligned data from Whole Genome Bisulphite Sequencing (WGBS) or Nanopore Sequencing datasets, performs a binomial test to help establish methylated sites, and formats the data for downstream differential methylation analyses. In addition to this, PryfynMeth can also: 

- Filter for coverage
- Reduce the number of sites to get tested in differential methylation analyses
- Combine all samples from the same treatment / tissue / sex to establish an overriding methylation pattern
- Generate PCAs

For further information, please examine the [Wiki documentation](https://github.com/C-L-Thomas/PryfynMeth/wiki/Home). The full pipeline for Nanopore sequencing can be found in the [Nanopore Section of the Wiki](https://github.com/C-L-Thomas/PryfynMeth/wiki/1.-Nanopore-Workflow).

# Installation
Download the repository:

```
git clone https://github.com/C-L-Thomas/PryfynMeth.git

```

Install:

```
cd PryfynMeth
pip install .
```

Test the installation works with the following command:

```
python PryfynMeth/pryfynmeth/binomial.py \
  -meta PryfynMeth/test_data/binomial/metadata.txt \
  -platform illu \
  -i PryfynMeth/test_data/binomial/input_dir \
  -o test
```



# Preprocessing
PryfynMeth can accept three input file types; Nanopore [MethylBed](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Nanopore_WorkFlow/example_methyl.bed), [Bismark bismark_methylation_extractor reports (bisulphite stranded)](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Bisulphite_Data/example_stranded.report.txt), and [Bismark coverage2cytosine cov files (bisulphite destranded)](https://github.com/C-L-Thomas/PryfynMeth/blob/main/Bisulphite_Data/example_destranded.cov). 

To find the full pipelines for each of these file types explore the [Wiki](https://github.com/C-L-Thomas/PryfynMeth/wiki/). The Nanopore wiki can be found [here](https://github.com/C-L-Thomas/PryfynMeth/wiki/1.-Nanopore-Workflow).



# Binomial Test 
To determine whether an individual site is methylated in whole genome methylation sequencing, it's common practice to assess if the observed methylation level is significantly higher than what would be expected by chance. This is typically done by comparing the observed methylation proportion to a statistical threshold, often using a binomial test. In whole genome bisulphite sequencing and Nanopore sequencing, this threshold is determined by the percentage methylation found in lambda spiked DNA. To identify methylated sites, you can use the `binomial.py` command. 



# Filtering
Once you have generated binomial results, you may wish to filter samples with low read count. The filter.py command does just that. By setting a threshold it will output three folders. The first will be a full list of your sites in your binomial test output, but with the values and statistics adjusted considering your filtering. The second will be a folde of each sample's methylated sites (sites with FDR less than 0.05). The final, will scan each of your methylated site files, and will add genomic locations that have been excluded from other samples, making sure each sample has the same number of input site.

```
python pryfynmeth/filter.py  -i Binomial_Results_Folder -f output_adjusted_sites -m output_pass_binomial -s output_shared_sites -threshold 10 -revert

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

# Compatibility with Differential Methylation Analyses

Currently, the only differential methylation toolkit that PryfynMeth is streamlined with is DSS. But if you have any additional requests please let me know.

# Errors & Requests

Whilst these scripts have been extensively trialed, errors may still occur. If you get any error messages, or any of the descriptions are unclear, please email **Christianluthomas@gmail.com**. Additionally, if you have any requests for additions to the pipeline feel free to email.

# Citations
