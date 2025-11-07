# Background
Welcome to the package PryfynMeth. PryfynMeth is a combination of tools designed to help the analysis of insect methylation data. It takes aligned data from Whole Genome Bisulphite Sequencing (WGBS) or Nanopore Sequencing datasets, performs a binomial test to help establish methylated sites, and formats the data for downstream differential methylation analyses. In addition to this, PryfynMeth can also: 

- Filter for coverage
- Reduce the number of sites to get tested in differential methylation analyses
- Combine all samples from the same treatment / tissue / sex to establish an overriding methylation pattern
- Generate [PCAs](https://github.com/C-L-Thomas/PryfynMeth/wiki/5.-PCA)
- Generate [Alignment Statistics](https://github.com/C-L-Thomas/PryfynMeth/wiki/6.-Statistics)
- Perform base letter [motif enrichment](https://github.com/C-L-Thomas/PryfynMeth/wiki/7.-Motif-Enrichment)

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
The preprocessing step allows users to determine the number of input sites for the binomial and consequentially determine the number of tests. It also formats input data ready for binomial testing.

# Binomial Test 
To determine whether an individual site is methylated in whole genome methylation sequencing, it's common practice to assess if the observed methylation level is significantly higher than what would be expected by chance. This is typically done by comparing the observed methylation proportion to a statistical threshold, often using a binomial test. In whole genome bisulphite sequencing and Nanopore sequencing, this threshold is determined by the percentage methylation found in lambda spiked DNA. To identify methylated sites, you can use the `binomial.py` command. 

# Filtering
Once you have generated binomial results, you may wish to filter samples with low read count. The filter.py command does just that. By setting a threshold it will output three folders. The first will be a full list of your sites in your binomial test output, but with the values and statistics adjusted considering your filtering. The second will be a folder of each sample's methylated sites (sites with FDR less than 0.05). The final, will scan each of your methylated site files, and will add genomic locations that have been excluded from other samples, making sure each sample has the same number of input site.

# Compatibility with Differential Methylation Analyses
Currently, the only differential methylation toolkit that PryfynMeth is streamlined with is DSS. But if you have any additional requests please let me know.

# Errors & Requests

Whilst these scripts have been extensively trialed, errors may still occur. If you get any error messages, or any of the descriptions are unclear, please email **Christianluthomas@gmail.com**. Additionally, if you have any requests for additions to the pipeline feel free to email.

# Citations
For citations please use:

Thomas, C., & Mallon, E. (2025). PryfynMeth: A pipeline for invertebrate methylation analysis [Report]. Zenodo. https://doi.org/10.5281/zenodo.17552905
