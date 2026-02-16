#!/usr/bin/env python3
"""
Calculate coverage statistics for chromosome-scale scaffolds.

This script calculates mean coverage and percentiles for each chromosome-scale
scaffold to establish the expected coverage distribution for the assembly.
"""

import argparse
import pysam
import sys
import numpy as np
import json


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate coverage statistics for chromosome scaffolds"
    )
    parser.add_argument(
        "--bam",
        required=True,
        help="Input BAM file with reads mapped to assembly"
    )
    parser.add_argument(
        "--chr_list",
        required=True,
        help="File with list of chromosome scaffold names (one per line)"
    )
    parser.add_argument(
        "--output_json",
        required=True,
        help="Output JSON file with coverage statistics"
    )
    parser.add_argument(
        "--output_tsv",
        required=True,
        help="Output TSV file with per-chromosome coverage"
    )
    
    return parser.parse_args()


def load_chr_list(chr_list_file):
    """Load list of chromosome names."""
    chr_names = []
    with open(chr_list_file) as f:
        for line in f:
            line = line.strip()
            if line:
                chr_names.append(line)
    return chr_names


def calculate_chr_coverage(bam_file, chr_names):
    """
    Calculate coverage statistics for chromosome scaffolds.
    
    Returns:
        dict: {chr_name: {'mean': X, 'median': Y, 'p25': Z, 'p75': W}}
        dict: {'genome_mean': X, 'genome_median': Y, 'genome_p25': Z, 'genome_p75': W}
    """
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    
    chr_stats = {}
    all_coverage_values = []
    
    for chr_name in chr_names:
        try:
            # Get chromosome length
            chr_length = bamfile.get_reference_length(chr_name)
            
            print(f"  Processing {chr_name} ({chr_length:,} bp)...", file=sys.stderr)
            
            # Get coverage array for entire chromosome
            # This returns tuple of 4 arrays (A, C, G, T counts)
            coverage_arrays = bamfile.count_coverage(
                chr_name,
                0,
                chr_length,
                quality_threshold=0
            )
            
            # Sum across all bases at each position
            coverage_per_base = np.sum([np.array(arr) for arr in coverage_arrays], axis=0)
            
            # Calculate statistics
            mean_cov = float(np.mean(coverage_per_base))
            median_cov = float(np.median(coverage_per_base))
            p25 = float(np.percentile(coverage_per_base, 25))
            p75 = float(np.percentile(coverage_per_base, 75))
            p10 = float(np.percentile(coverage_per_base, 10))
            p90 = float(np.percentile(coverage_per_base, 90))
            
            chr_stats[chr_name] = {
                'mean': mean_cov,
                'median': median_cov,
                'p10': p10,
                'p25': p25,
                'p75': p75,
                'p90': p90,
                'length': chr_length
            }
            
            # Add to genome-wide distribution
            all_coverage_values.extend(coverage_per_base.tolist())
            
            print(f"    Mean: {mean_cov:.1f}x, Median: {median_cov:.1f}x, "
                  f"25-75th: {p25:.1f}-{p75:.1f}x", file=sys.stderr)
            
        except Exception as e:
            print(f"  Warning: Could not process {chr_name}: {e}", file=sys.stderr)
            continue
    
    # Calculate genome-wide statistics
    if all_coverage_values:
        all_cov_array = np.array(all_coverage_values)
        genome_stats = {
            'mean': float(np.mean(all_cov_array)),
            'median': float(np.median(all_cov_array)),
            'p10': float(np.percentile(all_cov_array, 10)),
            'p25': float(np.percentile(all_cov_array, 25)),
            'p75': float(np.percentile(all_cov_array, 75)),
            'p90': float(np.percentile(all_cov_array, 90))
        }
    else:
        genome_stats = {}
    
    bamfile.close()
    
    return chr_stats, genome_stats


def write_outputs(chr_stats, genome_stats, output_json, output_tsv):
    """Write coverage statistics to output files."""
    
    # Write JSON with all stats
    output_data = {
        'chromosomes': chr_stats,
        'genome': genome_stats
    }
    
    with open(output_json, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    # Write TSV for easy viewing
    with open(output_tsv, 'w') as f:
        f.write("chromosome\tmean_coverage\tmedian_coverage\tp25\tp75\tp10\tp90\tlength\n")
        
        for chr_name in sorted(chr_stats.keys()):
            stats = chr_stats[chr_name]
            f.write(f"{chr_name}\t{stats['mean']:.2f}\t{stats['median']:.2f}\t"
                   f"{stats['p25']:.2f}\t{stats['p75']:.2f}\t{stats['p10']:.2f}\t"
                   f"{stats['p90']:.2f}\t{stats['length']}\n")
        
        # Add genome-wide row
        if genome_stats:
            f.write(f"GENOME_WIDE\t{genome_stats['mean']:.2f}\t{genome_stats['median']:.2f}\t"
                   f"{genome_stats['p25']:.2f}\t{genome_stats['p75']:.2f}\t"
                   f"{genome_stats['p10']:.2f}\t{genome_stats['p90']:.2f}\t-\n")


def main():
    """Main function."""
    args = parse_arguments()
    
    print("Loading chromosome list...", file=sys.stderr)
    chr_names = load_chr_list(args.chr_list)
    print(f"  Found {len(chr_names)} chromosomes", file=sys.stderr)
    
    print("Calculating coverage statistics...", file=sys.stderr)
    chr_stats, genome_stats = calculate_chr_coverage(args.bam, chr_names)
    
    print("\nGenome-wide coverage statistics:", file=sys.stderr)
    if genome_stats:
        print(f"  Mean: {genome_stats['mean']:.1f}x", file=sys.stderr)
        print(f"  Median: {genome_stats['median']:.1f}x", file=sys.stderr)
        print(f"  25-75th percentile: {genome_stats['p25']:.1f}-{genome_stats['p75']:.1f}x", 
              file=sys.stderr)
        print(f"  10-90th percentile: {genome_stats['p10']:.1f}-{genome_stats['p90']:.1f}x", 
              file=sys.stderr)
    
    print(f"\nWriting outputs...", file=sys.stderr)
    write_outputs(chr_stats, genome_stats, args.output_json, args.output_tsv)
    print(f"  JSON: {args.output_json}", file=sys.stderr)
    print(f"  TSV: {args.output_tsv}", file=sys.stderr)
    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
