#!/usr/bin/env python3
"""
Calculate binned coverage for contigs.

This script calculates coverage in fixed-size bins (e.g., 1kb) for each contig
to enable coverage track visualization in dotplots.
"""

import argparse
import pysam
import sys
import json
import numpy as np


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate binned coverage for contigs"
    )
    parser.add_argument(
        "--bam",
        required=True,
        help="Input BAM file with reads mapped to assembly"
    )
    parser.add_argument(
        "--contig_map",
        required=True,
        help="TSV file mapping contig names to scaffold regions"
    )
    parser.add_argument(
        "--output_json",
        required=True,
        help="Output JSON file with binned coverage data"
    )
    parser.add_argument(
        "--bin_size",
        type=int,
        default=1000,
        help="Bin size in bp (default: 1000)"
    )
    
    return parser.parse_args()


def read_contig_map(contig_map_file):
    """
    Read contig to scaffold mapping from TSV file.
    
    Returns:
        List of tuples: (contig_name, scaffold_name, start, end)
    """
    contigs = []
    with open(contig_map_file, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                contig_name = parts[0]
                scaffold_name = parts[1]
                scaffold_start = int(parts[4])
                scaffold_end = int(parts[5])
                contigs.append((contig_name, scaffold_name, scaffold_start, scaffold_end))
    return contigs


def calculate_binned_coverage(bam_file, contigs, bin_size):
    """
    Calculate coverage in fixed-size bins for each contig.
    
    Args:
        bam_file: Path to BAM file
        contigs: List of (contig_name, scaffold_name, start, end) tuples
        bin_size: Size of bins in bp
        
    Returns:
        Dictionary: {contig_name: {'bins': [cov1, cov2, ...], 'mean': X}}
    """
    coverage_data = {}
    
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    
    for idx, (contig_name, scaffold_name, start, end) in enumerate(contigs, 1):
        if idx % 100 == 0:
            print(f"  Processed {idx}/{len(contigs)} contigs...", file=sys.stderr)
        
        try:
            contig_length = end - start
            
            # Get coverage array for this region
            coverage_arrays = bamfile.count_coverage(
                scaffold_name,
                start,
                end,
                quality_threshold=0
            )
            
            # Sum across all bases (A, C, G, T)
            coverage_per_base = np.sum([np.array(arr) for arr in coverage_arrays], axis=0)
            
            # Calculate binned coverage
            num_bins = (contig_length + bin_size - 1) // bin_size  # Ceiling division
            binned_coverage = []
            
            for bin_idx in range(num_bins):
                bin_start = bin_idx * bin_size
                bin_end = min((bin_idx + 1) * bin_size, contig_length)
                
                # Average coverage in this bin
                if bin_end > bin_start:
                    bin_cov = float(np.mean(coverage_per_base[bin_start:bin_end]))
                else:
                    bin_cov = 0.0
                
                binned_coverage.append(round(bin_cov, 2))
            
            # Calculate mean coverage for entire contig
            mean_cov = float(np.mean(coverage_per_base))
            
            coverage_data[contig_name] = {
                'bins': binned_coverage,
                'mean': round(mean_cov, 2),
                'bin_size': bin_size,
                'length': contig_length
            }
            
        except Exception as e:
            print(f"  Warning: Could not process {contig_name}: {e}", file=sys.stderr)
            continue
    
    bamfile.close()
    
    return coverage_data


def write_output(coverage_data, output_file):
    """Write binned coverage data to JSON file."""
    with open(output_file, 'w') as f:
        json.dump(coverage_data, f, indent=2)


def main():
    """Main function."""
    args = parse_arguments()
    
    print(f"Reading contig map from {args.contig_map}...", file=sys.stderr)
    contigs = read_contig_map(args.contig_map)
    print(f"  Found {len(contigs)} contigs", file=sys.stderr)
    
    print(f"Calculating binned coverage (bin size: {args.bin_size} bp)...", file=sys.stderr)
    coverage_data = calculate_binned_coverage(args.bam, contigs, args.bin_size)
    print(f"  Calculated coverage for {len(coverage_data)} contigs", file=sys.stderr)
    
    print(f"Writing output to {args.output_json}...", file=sys.stderr)
    write_output(coverage_data, args.output_json)
    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
