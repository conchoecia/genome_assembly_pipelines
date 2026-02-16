#!/usr/bin/env python3
"""
Calculate per-contig coverage from BAM file using contig regions.

This script calculates mean coverage for each contig by extracting coverage
from specific regions of scaffolds based on the contig-to-scaffold mapping.
"""

import argparse
import pysam
import sys
from collections import defaultdict


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate per-contig coverage from BAM file"
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
        "--output",
        required=True,
        help="Output TSV file with per-contig coverage"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads (not used currently)"
    )
    
    return parser.parse_args()


def read_contig_map(contig_map_file):
    """
    Read contig to scaffold mapping from TSV file.
    
    Args:
        contig_map_file: Path to TSV file with contig mapping
        
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


def calculate_coverage(bam_file, contigs):
    """
    Calculate mean coverage for each contig.
    
    Args:
        bam_file: Path to BAM file
        contigs: List of (contig_name, scaffold_name, start, end) tuples
        
    Returns:
        Dictionary mapping contig_name to mean_coverage
    """
    coverage_dict = {}
    
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    
    for contig_name, scaffold_name, start, end in contigs:
        try:
            # Get coverage array for this region
            # pysam uses 0-based coordinates
            coverage_array = bamfile.count_coverage(
                scaffold_name, 
                start, 
                end,
                quality_threshold=0
            )
            
            # Sum coverage across all bases (A, C, G, T)
            total_coverage = 0
            num_positions = end - start
            
            for pos_idx in range(num_positions):
                pos_coverage = sum(arr[pos_idx] for arr in coverage_array)
                total_coverage += pos_coverage
            
            # Calculate mean coverage
            if num_positions > 0:
                mean_coverage = total_coverage / num_positions
            else:
                mean_coverage = 0.0
            
            coverage_dict[contig_name] = mean_coverage
            
        except Exception as e:
            print(f"Warning: Could not calculate coverage for contig {contig_name}: {e}", 
                  file=sys.stderr)
            coverage_dict[contig_name] = 0.0
    
    bamfile.close()
    
    return coverage_dict


def write_output(coverage_dict, output_file):
    """
    Write coverage results to TSV file.
    
    Args:
        coverage_dict: Dictionary mapping contig_name to mean_coverage
        output_file: Path to output TSV file
    """
    with open(output_file, 'w') as f:
        f.write("contig_name\tmean_coverage\n")
        for contig_name in sorted(coverage_dict.keys()):
            mean_coverage = coverage_dict[contig_name]
            f.write(f"{contig_name}\t{mean_coverage:.2f}\n")


def main():
    """Main function."""
    args = parse_arguments()
    
    print(f"Reading contig map from {args.contig_map}...", file=sys.stderr)
    contigs = read_contig_map(args.contig_map)
    print(f"  Found {len(contigs)} contigs", file=sys.stderr)
    
    print(f"Calculating coverage from {args.bam}...", file=sys.stderr)
    coverage_dict = calculate_coverage(args.bam, contigs)
    print(f"  Calculated coverage for {len(coverage_dict)} contigs", file=sys.stderr)
    
    print(f"Writing output to {args.output}...", file=sys.stderr)
    write_output(coverage_dict, args.output)
    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
