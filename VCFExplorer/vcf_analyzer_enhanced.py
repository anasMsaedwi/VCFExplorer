import pysam
import numpy as np
import pandas as pd
from collections import defaultdict

def analyze_vcf(vcf_file):
    # Open the VCF file
    vcf = pysam.VariantFile(vcf_file)

    # Initialize counters and data structures
    total_variants = 0
    variant_types = {"SNP": 0, "Insertion": 0, "Deletion": 0, "Multi-allelic": 0}
    chromosome_counts = defaultdict(int)
    allele_frequencies = []
    genotype_counts = defaultdict(int)
    sample_heterozygosity = defaultdict(int)
    sample_missing_data = defaultdict(int)
    quality_scores = []
    depths = []

    # Iterate through each variant in the VCF file
    for record in vcf:
        total_variants += 1

        # Count variant types
        if len(record.ref) == 1 and all(len(alt) == 1 for alt in record.alts):
            variant_types["SNP"] += 1
        elif any(len(alt) > len(record.ref) for alt in record.alts):
            variant_types["Insertion"] += 1
        elif any(len(alt) < len(record.ref) for alt in record.alts):
            variant_types["Deletion"] += 1
        else:
            variant_types["Multi-allelic"] += 1

        # Count variants per chromosome
        chromosome_counts[record.chrom] += 1

        # Calculate allele frequency (AF)
        af = record.info.get("AF", [None])[0]
        if af is not None:
            allele_frequencies.append(af)

        # Collect quality scores and depth
        quality_scores.append(record.qual)
        dp = record.info.get("DP", None)
        if dp is not None:
            depths.append(dp)

        # Count genotypes and calculate sample-level statistics
        for sample_name, sample in record.samples.items():
            gt = sample.get("GT")
            if gt is not None:
                gt_str = "/".join(map(str, gt))
                genotype_counts[gt_str] += 1
                if gt_str == "0/1":
                    sample_heterozygosity[sample_name] += 1
            else:
                sample_missing_data[sample_name] += 1

    # Calculate statistics
    stats = {
        "total_variants": total_variants,
        "variant_types": variant_types,
        "chromosome_counts": dict(chromosome_counts),
        "mean_allele_frequency": np.mean(allele_frequencies) if allele_frequencies else None,
        "genotype_counts": dict(genotype_counts),
        "mean_quality_score": np.mean(quality_scores) if quality_scores else None,
        "mean_depth": np.mean(depths) if depths else None,
        "sample_heterozygosity": {sample: count / total_variants for sample, count in sample_heterozygosity.items()},
        "sample_missing_data": {sample: count / total_variants for sample, count in sample_missing_data.items()}
    }

    return stats

def write_results_to_file(stats, output_file):
    """Write the statistics to a text file."""
    with open(output_file, "w") as f:
        f.write("VCF Statistics Summary:\n")
        f.write("======================\n")
        f.write(f"Total Variants: {stats['total_variants']}\n")
        f.write("\nVariant Types:\n")
        for vt, count in stats["variant_types"].items():
            f.write(f"  {vt}: {count}\n")
        f.write("\nChromosome Distribution:\n")
        for chrom, count in stats["chromosome_counts"].items():
            f.write(f"  {chrom}: {count}\n")
        f.write(f"\nMean Allele Frequency: {stats['mean_allele_frequency']}\n")
        f.write(f"Mean Quality Score: {stats['mean_quality_score']}\n")
        f.write(f"Mean Depth: {stats['mean_depth']}\n")
        f.write("\nGenotype Counts:\n")
        for gt, count in stats["genotype_counts"].items():
            f.write(f"  {gt}: {count}\n")
        f.write("\nSample Heterozygosity Rates:\n")
        for sample, rate in stats["sample_heterozygosity"].items():
            f.write(f"  {sample}: {rate:.4f}\n")
        f.write("\nSample Missing Data Rates:\n")
        for sample, rate in stats["sample_missing_data"].items():
            f.write(f"  {sample}: {rate:.4f}\n")

def write_results_to_csv(stats, output_file):
    """Write the statistics to a CSV file."""
    with open(output_file, "w") as f:
        f.write("Statistic,Value\n")
        f.write(f"Total Variants,{stats['total_variants']}\n")
        for vt, count in stats["variant_types"].items():
            f.write(f"Variant Type - {vt},{count}\n")
        for chrom, count in stats["chromosome_counts"].items():
            f.write(f"Chromosome - {chrom},{count}\n")
        f.write(f"Mean Allele Frequency,{stats['mean_allele_frequency']}\n")
        f.write(f"Mean Quality Score,{stats['mean_quality_score']}\n")
        f.write(f"Mean Depth,{stats['mean_depth']}\n")
        for gt, count in stats["genotype_counts"].items():
            f.write(f"Genotype - {gt},{count}\n")
        for sample, rate in stats["sample_heterozygosity"].items():
            f.write(f"Sample Heterozygosity - {sample},{rate:.4f}\n")
        for sample, rate in stats["sample_missing_data"].items():
            f.write(f"Sample Missing Data - {sample},{rate:.4f}\n")

def main():
    # Path to the VCF file
    vcf_file = "/home/anas-saedwi/Documents/Python_progict/P_vcf/Sample_data.vcf"  # Replace with your VCF file path
    output_txt_file = "/home/anas-saedwi/Documents/Python_progict/P_vcf/vcf_statistics_summary.txt"  # Text output file
    output_csv_file = "/home/anas-saedwi/Documents/Python_progict/P_vcf/vcf_statistics_summary.csv"  # CSV output file

    # Analyze the VCF file
    stats = analyze_vcf(vcf_file)

    # Write results to text and CSV files
    write_results_to_file(stats, output_txt_file)
    write_results_to_csv(stats, output_csv_file)

    # Print confirmation
    print(f"Results written to {output_txt_file} and {output_csv_file}")

if __name__ == "__main__":
    main()
