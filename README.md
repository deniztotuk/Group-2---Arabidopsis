# Group-2---Arabidopsis

#### Analysis of Polyploidy in Arabidopsis Species

## Overview
This project explores the phenomenon of polyploidy in two Arabidopsis species: *Arabidopsis arenosa* and *Arabidopsis lyrata*. Our research focuses on understanding how independently formed autopolyploids in these species may have exchanged adaptive alleles , potentially influencing their evolutionary trajectories.

### Objectives
- To analyze allele frequency spectra across normal *arenosa* populations, autohexaploid *lyrata* populations, and their hybrids to determine differences and implications for survival and diversity.
- To explore the possibility of hybrids exhibiting greater genetic diversity, comparing them to pure populations and assessing if they resemble pseudo-allopolyploids.

### Installation
This project requires Python, R and several dependencies. Many of the dependencies used in this project can be easily installed via pip command for all common operating systems by using terminal or PowerShell.

If you are working on windows machine before using pip command installing python3.12.1 essential and can be easily installed from: [(https://www.python.org/downloads/release/python-3121/)](https://www.python.org/downloads/release/python-3121/)

In this project we have used R statistical programming language. R version 4.3.2 can be installed from: [(https://cran.r-project.org/bin/windows/base/old/4.3.2/)](https://cran.r-project.org/bin/windows/base/old/4.3.2/)

```bash
# Install required Python packages
pip install matplotlib==3.8.3
pip install numpy==1.26.0
pip install dadi==2.3.3
pip install pysam==0.22.0
pip install pandas==2.2.2
```
Ins
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
Python script for converting VCF to Phylip format for doing phylogenetic analysis using SplitsTree. Script is created by Deniz Totuk.
```python
# Import necessary libraries: pysam for handling VCF files and sys for accessing command-line arguments.
import pysam
import sys

# Define the function vcf_to_phylip which takes two parameters: the input VCF file and the output Phylip file.
def vcf_to_phylip(vcf_file, output_file):
    # Open the VCF file using pysam.
    vcf = pysam.VariantFile(vcf_file)
    # Extract sample names from the VCF file's header.
    samples = list(vcf.header.samples)
    # Initialize a dictionary to hold sequence data for each sample.
    sequence_data = {sample: [] for sample in samples}

    # Iterate over each record (variant) in the VCF file.
    for record in vcf:
        # Store the reference and all alternate alleles of the variant.
        alleles = [record.ref] + list(record.alts)
        # Iterate over each sample to determine the genotype and corresponding allele sequence.
        for sample in samples:
            genotype = record.samples[sample]['GT']
            # Check if genotype information is missing and handle it by adding 'N' (common symbol for unknown).
            if None in genotype:  # Handle missing data
                sequence_data[sample].append('N')
            else:
                # Add the allele corresponding to the first allele of the genotype tuple.
                sequence_data[sample].append(alleles[genotype[0]])

    # Open the output file for writing.
    with open(output_file, 'w') as f:
        # Write the header of the Phylip file, which includes the number of samples and the length of the sequences.
        f.write(f"{len(samples)} {len(sequence_data[samples[0]])}\n")
        # Write the sequence data for each sample.
        for sample, sequence in sequence_data.items():
            f.write(f"{sample} {''.join(sequence)}\n")

# Check if the script is executed directly (i.e., not imported as a module), and if so, call the function with command-line arguments.
if __name__ == "__main__":
    vcf_to_phylip(sys.argv[1], sys.argv[2])  # sys.argv[1] is the input VCF file, sys.argv[2] is the output Phylip file
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
This script describes how to run the SplitsTree software for phylogenetic analysis using a shell script. Script is created by Kavithi Jayasundara.
```bash
#!/bin/bash
# Command to run SplitsTree
splitstree -g -i $final_output.phyphilp -o $final_output.nexus
```


# Site Frequency Spectrum Calculation
This script calculates and plots the Site Frequency Spectrum (SFS) for each species using the dadi Python package. Script is created by Deniz Totuk.

```python
# Import necessary libraries for numerical operations, demographic analysis, and plotting
import numpy as np
import dadi
import matplotlib.pyplot as plt

def read_allele_counts(file_path):
    """Reads allele counts from a file; each line is an allele count at a site."""
    # Open the file for reading.
    with open(file_path, 'r') as f:
        # Read lines from the file, strip whitespace, and convert each line to an integer if not empty
        counts = [int(line.strip()) for line in f if line.strip()]
    return counts  # Return the list of allele counts

def calculate_sfs(allele_counts, projection):
    """Calculates the unfolded Site Frequency Spectrum (SFS) from allele counts."""
    # Convert the list of allele counts into a dadi Spectrum object
    fs = dadi.Spectrum(allele_counts)
    # Project the spectrum to a smaller sample size to smooth the spectrum
    fs_projected = fs.project([projection])
    return fs_projected  # Return the projected spectrum

def plot_sfs(sfs, label):
    """Plots the SFS."""
    # Plot the SFS using a line plot with markers, skipping the zero index as it's often not informative
    plt.plot(range(1, len(sfs)), sfs[1:], label=label, marker='o')

# Read allele counts from files for two populations
arenosa_counts = read_allele_counts('arenosa_for_sfs.txt')
lyrata_counts = read_allele_counts('lyrata_for_sfs.txt')

# Define the size to which the allele counts are projected
projection_size = 20
# Calculate the SFS for each population using the defined projection size
arenosa_sfs = calculate_sfs(arenosa_counts, projection_size)
lyrata_sfs = calculate_sfs(lyrata_counts, projection_size)

# Set up the plotting area with a specific size
plt.figure(figsize=(8, 6))
# Plot the SFS for each population
plot_sfs(arenosa_sfs, 'Arenosa')
plot_sfs(lyrata_sfs, 'Lyrata')
# Set the labels for the axes and the plot title
plt.xlabel('Derived allele frequency')
plt.ylabel('Count')
plt.title('Site Frequency Spectrum (SFS)')
# Show a legend to identify the populations
plt.legend()
# Adjust plot layout for better spacing
plt.tight_layout()
# Display the plot
plt.show()
```

## References

1. Harris, Charles R., et al. "Array programming with NumPy." Nature 585.7825 (2020): 357-362. Available at: [https://numpy.org](https://numpy.org)
2. McKinney, Wes. "Data Structures for Statistical Computing in Python." Proceedings of the 9th Python in Science Conference. 2010. Available at: [https://pandas.pydata.org](https://pandas.pydata.org)
3. Hunter, John D. "Matplotlib: A 2D Graphics Environment." Computing in Science & Engineering 9.03 (2007): 90-95. Available at: [https://matplotlib.org](https://matplotlib.org)
5. Gutenkunst, Ryan N., et al. "Inferring the joint demographic history of multiple populations from multidimensional SNP frequency data." PLoS Genetics 5.10 (2009): e1000695. Available at: [https://dadi.readthedocs.io/en/latest/](https://dadi.readthedocs.io/en/latest/)
6. Heger, Andreas, and Sean R. Eddy. "pysam: a Python library for working with SAM/BAM/CRAM formats." Available at: [https://github.com/pysam-developers/pysam](https://github.com/pysam-developers/pysam)
7. Van Rossum, Guido, and Drake, Fred L. Jr., "Python 3 Reference Manual", Python Software Foundation, 2009. Available at: [https://docs.python.org/3/library/sys.html](https://docs.python.org/3/library/sys.html)
8. R Core Team. "R: A language and environment for statistical computing." R Foundation for Statistical Computing, Vienna, Austria, 2021. Available at: [https://www.R-project.org/](https://www.R-project.org/)
9. RStudio Team. "RStudio: Integrated Development Environment for R." RStudio, PBC, Boston, MA, 2021. Available at: [https://www.rstudio.com/](https://www.rstudio.com/)
10. Li, Heng, et al. "The Sequence Alignment/Map format and SAMtools." Bioinformatics 25.16 (2009): 2078-2079. Available at: [https://academic.oup.com/bioinformatics/article/25/16/2078/204688](https://academic.oup.com/bioinformatics/article/25/16/2078/204688)
- "bcftools: Utilities for variant calling and manipulating VCFs and BCFs." Available at: [http://www.htslib.org/doc/bcftools.html](http://www.htslib.org/doc/bcftools.html)

