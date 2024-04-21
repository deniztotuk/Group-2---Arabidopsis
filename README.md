# Group-2---Arabidopsis

#### Analysis of Polyploidy in Arabidopsis Species

## Overview
This project explores the phenomenon of polyploidy in two Arabidopsis species: *Arabidopsis arenosa* and *Arabidopsis lyrata*. Our research focuses on understanding how independently formed autopolyploids in these species may have exchanged adaptive alleles , potentially influencing their evolutionary trajectories.

### Objectives
- To analyze allele frequency spectra across normal *arenosa* populations, autohexaploid *lyrata* populations, and their hybrids to determine differences and implications for survival and diversity.
- To explore the possibility of hybrids exhibiting greater genetic diversity, comparing them to pure populations and assessing if they resemble pseudo-allopolyploids.

## Installation
This project requires Python, R and several dependencies:

```bash
# Install required Python packages
pip install matplotlib==3.4.3
pip install numpy==1.19.5
# Add additional dependencies as necessary
```
### Data
This project uses VCF files from 632 Arabidopsis arenosa individuals and 272 Arabidopsis lyrata individuals (including hybrids). The data is analyzed to understand genetic diversity and differentiation metrics.

# Data Formats and Requirements:
For input and output we primarily handle VCF files. Please ensure all filenames follow the specific conventions as outlined in the scripts.

### Scripts and Usage:

# Initial Data Preparation
Filtering VCF Files
The initial VCF files are filtered using GATK's SelectVariants to include only specific samples. This step is crucial for reducing data complexity and focusing analysis on relevant populations. Script is created by 

```bash
# Activate the GATK environment
conda activate /shared/conda/gatk

# Run GATK SelectVariants to filter VCF for specified samples
gatk SelectVariants \
-V /workhere/students_2023/Levi_resources/Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf \
-sn KAG-01tl -sn KAG-02tl -sn KAG-03tl -sn KAG-04tl -sn KAG-05tl \
-sn MOD-01tl -sn MOD-02tl -sn MOD-03tl -sn MOD-04tl \
-sn LIC-01tl -sn LIC-02tl -sn LIC-03tl \
-sn FRE-01tl -sn FRE-02tl -sn FRE-03tl -sn FRE-04tl -sn FRE-07tl \
-sn OCH-01tl -sn OCH-02tl -sn OCH-03tl -sn OCH-04tl -sn OCH-06tl -sn OCH-07tl -sn OCH-08tl \
-sn HAB-01tl -sn HAB-02tl -sn HAB-03tl \
-sn KEH-01tl -sn KEH-04tl -sn KEH-05tl -sn KEH-10tl \
-O /workhere/students_2023/Levi_resources/group2/filtered_2.vcf
```


