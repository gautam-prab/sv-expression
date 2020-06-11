# SV Expression Readme

This repository is a collection of scripts used to analyze the effects of structural variants on gene expression using data from the 1000 Genomes Project. All of the scripts have a `--help` option for a list of command-line options.

### Table of Contents
- [Summary Statistics for a VCF](#vcf)
- [Finding Significant SV-Gene Associations](#svgene)
- [Plots with SV-Gene Associations](#geneplots)
- [Permutation Tests for Enrichment](#permutation)
- [Packages Used](#packages)

<a name="vcf"></a>
## Summary Statistics for a Structural Variants

Given a VCF with per-sample genotypes for a set of structural variants, `process_vcf.py` allows the generation of data that can be plotted for some basic summary statistics.

#### Pre-Processing the VCF

The script `preprocess_vcf.py` filters SVs based on HWE (Hardy-Weinberg Equilibrium p-value), # uncalled, and minor allele frequency:

```
python pre_process_vcf.py file.vcf.gz --outfile filtered.vcf.gz
```

the help options `--uncalled`, `--hwe`, `--maf`, and `--ref` allow for specific filtering changes.

#### Statistics Plots

There are two possible output files: lengthfile, which is a python dictionary containing the lengths of all SVs for each SV type, and samplefile, which is a CSV of each sample with the number of each type of SV that they have.

```
python process_vcf.py --vcffile input.vcf.gz --lengthfile length.dict --samplefile samples.csv
```

These can be plotted:

```
python make_svsize_plot.py --dictfile length.dict
```

for a plot of SV sizes, and

```
python make_svpopulation_plot.py --csvfile samples.csv --population super_pop
```

where the population argument separates samples by a particular population (which should be specified in the --populationfile argument).

<a name="svgene"></a>
## Finding Significant SV-Gene Associations
Using a VCF and an RNA-seq file, I find significant SV-gene associations using `process_rnaseq.py` and `process_pairs.py`.

#### Data Specifications

The VCF should be bgzipped and indexed:
`bgzip file.vcf && tabix -p vcf file.vcf.gz`.

The RNA-seq file should be a TSV arranged as follows:
```
TargetID    Gene_Symbol Chr Coord   Sample1 Sample2 ...
ENSG0001    ENSG0001    1   100     0.412   1.323
ENSG0002    ENSG0002    1   40      0.451   2.31
...
```
These data should be PEER-normalized and provided in RPKM.

#### Association Tests and Effect Size Calculation

Use `process_rnaseq.py` to do association tests:

```
python process_rnaseq.py --vcffile file.vcf.gz --rnafile file.txt --p X --pfile out_pfile.csv --outfile out.vcf.gz
```

If you need to decide the correct p-value threshold, you can use `python plot_pvals.py out_pfile.csv` to get appropriate thresholds for Bonferroni or false discovery rate (FDR) correction. Then run `process_rnaseq.py` again using this new p-value threshold.

This associates variants with a gene field in the INFO column of `out.vcf.gz`, then use `process_pairs.py` on this file to calculate beta, standard error, and R<sup>2</sup>:

```
python process_pairs.py file.vcf.gz --outcsv out.csv --outvcf out.vcf
```

You can get the output here in either VCF (within INFO column) or CSV format.

<a name="geneplots"></a>
## Plots with SV-Gene Associations

#### Box Plots of Beta Values
To do
`python make_beta_boxplots.py`

#### Violin Plot of Beta Values
To do
`python make_violin_plots.py`

<a name="permutation"></a>
## Permutation Tests for Enrichments
To do

<a name="packages"></a>
## Packages Used
- `numpy`, `pandas`, `scipy`, `sklearn`, `matplotlib`
- [`cyvcf2`](https://github.com/brentp/cyvcf2)
- [`snphwe`](https://github.com/jeremymcrae/snphwe)
- [`pybedtools`](https://github.com/daler/pybedtools)
- [`pyliftover`](https://github.com/konstantint/pyliftover)
- [`mygene`](https://github.com/biothings/mygene.info)
- [`intervaltree`](https://github.com/chaimleib/intervaltree)
