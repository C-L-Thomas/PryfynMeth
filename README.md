# Background
Welcome to the package PryfynMeth. PryfynMeth is a combination of tools designed to help the analysis of insect methylation data. It takes aligned data from Whole Genome Bisulphite Sequencing (WGBS) or Nanopore Sequencing datasets, performs a binomial test to help establish methylated sites, and formats the data for downstream differential methylation analyses. In addition to this, PryfynMeth can also: 

- Filter for coverage
- Reduce the number of sites to get tested in differential methylation analyses
- Combine all samples from the same treatment / tissue / sex to establish an overriding methylation pattern
- Generate PCAs

For information on the full use of each script, please examine the [Wiki documentation](https://github.com/C-L-Thomas/PryfynMeth/wiki/Home). The full pipeline for Nanopore sequencing can be found in the [Nanopore Section of the Wiki](https://github.com/C-L-Thomas/PryfynMeth/wiki/1.-Nanopore-Workflow). The full pipeline for WGBS Stranded data can be found in the [Stranded WGBS Workflow Section of the Wiki](https://github.com/C-L-Thomas/PryfynMeth/wiki/1.-Stranded-WGBS-Workflow). The full pipeline for WGBS Destranded data can be found in the [Destranded WGBS Workflow Section of the Wiki](https://github.com/C-L-Thomas/PryfynMeth/wiki/3.-Detranded-WGBS-Workflow).

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


# Compatibility with Differential Methylation Analyses

Currently, the only differential methylation toolkit that PryfynMeth is streamlined with is DSS. But if you have any additional requests please let me know.

# Errors & Requests

Whilst these scripts have been extensively trialed, errors may still occur. If you get any error messages, or any of the descriptions are unclear, please email **Christianluthomas@gmail.com**. Additionally, if you have any requests for additions to the pipeline feel free to email.

# Citations
