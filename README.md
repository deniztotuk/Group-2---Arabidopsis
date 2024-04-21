# Group-2---Arabidopsis

#### Analysis of Polyploidy in Arabidopsis Species

## Overview
This project explores the phenomenon of polyploidy in two Arabidopsis species: *Arabidopsis arenosa* and *Arabidopsis lyrata*. Our research focuses on understanding how independently formed autopolyploids in these species may have exchanged adaptive alleles , potentially influencing their evolutionary trajectories.

### Objectives
- To analyze allele frequency spectra across normal *arenosa* populations, autohexaploid *lyrata* populations, and their hybrids to determine differences and implications for survival and diversity.
- To explore the possibility of hybrids exhibiting greater genetic diversity, comparing them to pure populations and assessing if they resemble pseudo-allopolyploids.

### Installation
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
The initial VCF files are filtered using GATK's SelectVariants to include only specific samples. This step is crucial for reducing data complexity and focusing analysis on relevant populations. Script is created by Jasmin Kuar.

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
Identify common genomic sites between two datasets. Script is created by Kavithi Jayasundara.
```bash
Copy code
awk 'NR==FNR {common[$1,$2]; next} ($1,$2) in common' lyrata_272_with_some_hybrids.txt arenosa_632.txt > common_sites.txt
```
Python script for converting VCF to Phylip format for doing phylogenetic analysis using SplitsTree.
```python
import pysam
import sys

def vcf_to_phylip(vcf_file, output_file):
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    sequence_data = {sample: [] for sample in samples}

    for record in vcf:
        alleles = [record.ref] + list(record.alts)
        for sample in samples:
            genotype = record.samples[sample]['GT']
            if None in genotype:  # Handle missing data
                sequence_data[sample].append('N')
            else:
                sequence_data[sample].append(alleles[genotype[0]])

    with open(output_file, 'w') as f:
        f.write(f"{len(samples)} {len(sequence_data[samples[0]])}\n")
        for sample, sequence in sequence_data.items():
            f.write(f"{sample} {''.join(sequence)}\n")

if __name__ == "__main__":
    vcf_to_phylip(sys.argv[1], sys.argv[2])
```

Command for converting VCF to Phylip format for doing phylogenetic analysis using SplitsTree.
```bash
chmod +x vcf_to_phylip.py
./vcf2phylip.py -i fully_final_filtered.vcf.gz -n --output-prefix final_output
```

# Allele Frequency Histogram
This script calculates and plots the allele frequency histogram for the BZD population. Ensure you have the fully_final_filtered.vcf.gz ready. Script is created by Kavithi Jayasundara.

```bash
# Extract individual and population names
bcftools query -l fully_final_filtered.vcf.gz | awk -F '-' '{print $0 "\t" $1}' | sort -k1,1 > pops_complete.txt

# Compile and run the frequency analysis
gcc poly_freq.c -o poly_freq -lm
./poly_freq -vcf fully_final_filtered.vcf -pops pops_complete.txt > info_final.tsv

# Plot the histogram using R
df <- read.table(file ='info_final.tsv', header = TRUE, sep = '\t')
ggplot(data = df, aes(x = BZD)) +
  geom_histogram(color = 'black', fill = 'white', bins = 10, breaks = seq(0, 1, by = 0.1)) +
  labs(title = "Allele Frequency Histogram for BZD Population", x = "Allele Frequency", y = "F
```
# Phylogenetic Analysis Using Splits Tree
This script describes how to run the SplitsTree software for phylogenetic analysis using a shell script.
```bash
#!/bin/bash
# Command to run SplitsTree
splitstree -g -i $INPUT_PHYLIP -o $OUTPUT_NEXUS
```


# Site Frequency Spectrum Calculation
This script calculates and plots the Site Frequency Spectrum (SFS) for each species using the dadi Python package. Script is created by Deniz Totuk.

```python
import numpy as np
import dadi
import matplotlib.pyplot as plt

def read_allele_counts(file_path):
    """Reads allele counts from a file; each line is an allele count at a site."""
    with open(file_path, 'r') as f:
        counts = [int(line.strip()) for line in f if line.strip()]
    return counts

def calculate_sfs(allele_counts, projection):
    """Calculates the unfolded SFS from allele counts."""
    fs = dadi.Spectrum(allele_counts)
    fs_projected = fs.project([projection])
    return fs_projected

def plot_sfs(sfs, label):
    """Plots the SFS."""
    plt.plot(range(1, len(sfs)), sfs[1:], label=label, marker='o')

# Example usage
arenosa_counts = read_allele_counts('arenosa_for_sfs.txt')
lyrata_counts = read_allele_counts('lyrata_for_sfs.txt')
projection_size = 20
arenosa_sfs = calculate_sfs(arenosa_counts, projection_size)
lyrata_sfs = calculate_sfs(lyrata_counts, projection_size)

plt.figure(figsize=(8, 6))
plot_sfs(arenosa_sfs, 'Arenosa')
plot_sfs(lyrata_sfs, 'Lyrata')
plt.xlabel('Derived allele frequency')
plt.ylabel('Count')
plt.title('Site Frequency Spectrum (SFS)')
plt.legend()
plt.tight_layout()
plt.show()
```



