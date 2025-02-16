# VCFExplorer
VCFAnalyzer VCFAnalyzer is a Python-based tool for analyzing and processing Variant Call Format (VCF) files. It provides detailed statistics on variant types (SNPs, insertions, deletions, and multi-allelic variants), allele frequencies, genotype distributions, chromosome-wise variant counts, quality scores, and read depths. 



# VCF File Analyzer

This Python script analyzes VCF (Variant Call Format) files and generates statistics including:
- Total variants
- Variant types (SNPs, Insertions, Deletions, Multi-allelic)
- Chromosome distribution
- Allele frequencies
- Quality scores
- Genotype counts
- Sample heterozygosity and missing data rates

## Features
- Detailed statistics summary saved as a text file
- CSV export for easy data analysis
- Handles large VCF files efficiently

## Requirements
- Python 3.x
- Dependencies: `pysam`, `numpy`, `pandas`

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/your-username/P_vcf.git
    cd P_vcf
    ```

2. Install the required packages:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

Place your VCF file in the `data/` directory and specify the file path in `analyze_vcf.py`. Then, run:

```bash
python analyze_vcf.py

