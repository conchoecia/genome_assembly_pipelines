#!/usr/bin/env python3
"""
Parse scaffold-to-scaffold mapping BAM file to calculate mapping statistics.

This script analyzes how much each non-chromosome scaffold maps to each 
chromosome-scale scaffold. It outputs a TSV file with mapping percentages
for each scaffold to each chromosome.
"""

import argparse
import pysam
from Bio import SeqIO
from collections import defaultdict
import sys


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse scaffold mapping BAM to calculate mapping statistics"
    )
    parser.add_argument(
        "--bam",
        required=True,
        help="Input BAM file (non-chr scaffolds mapped to chr scaffolds)"
    )
    parser.add_argument(
        "--non_chr_fasta",
        required=True,
        help="FASTA file containing non-chromosome scaffolds"
    )
    parser.add_argument(
        "--chr_list",
        required=True,
        help="Text file with list of chromosome scaffold names (one per line)"
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=0,
        help="Minimum mapping quality to consider (default: 0)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file with mapping statistics"
    )
    
    return parser.parse_args()


def get_scaffold_lengths(fasta_file):
    """
    Get lengths of all scaffolds from FASTA file.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        Dictionary mapping scaffold name to length
    """
    lengths = {}
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            lengths[record.id] = len(record.seq)
    return lengths


def get_chromosome_list(chr_list_file):
    """
    Read chromosome list from file.
    
    Args:
        chr_list_file: Path to text file with chromosome names
        
    Returns:
        List of chromosome names
    """
    chromosomes = []
    with open(chr_list_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                chromosomes.append(line)
    return chromosomes


def calculate_mapping_stats(bam_file, scaffold_lengths, chromosomes, min_mapq):
    """
    Calculate mapping statistics from BAM file.
    
    For each non-chr scaffold, calculate how many bases map to each chr scaffold.
    
    Args:
        bam_file: Path to BAM file
        scaffold_lengths: Dictionary of scaffold name to length
        chromosomes: List of chromosome names
        min_mapq: Minimum mapping quality to consider
        
    Returns:
        Dictionary with mapping statistics for each scaffold
    """
    # Initialize data structures
    # mapping_stats[scaffold_name][chr_name] = aligned_bases
    mapping_stats = defaultdict(lambda: defaultdict(int))
    
    # Process BAM file
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    
    for read in bamfile.fetch():
        # Skip unmapped reads
        if read.is_unmapped:
            continue
            
        # Skip reads below mapping quality threshold
        if read.mapping_quality < min_mapq:
            continue
        
        query_name = read.query_name  # non-chr scaffold name
        ref_name = read.reference_name  # chr scaffold name
        
        # Calculate aligned bases (using alignment blocks to skip insertions/deletions)
        # get_blocks() returns list of (start, end) tuples
        aligned_bases = sum(end - start for start, end in read.get_blocks())
        
        mapping_stats[query_name][ref_name] += aligned_bases
    
    bamfile.close()
    
    # Calculate percentages and statistics
    results = []
    
    for scaffold_name in sorted(scaffold_lengths.keys()):
        scaffold_length = scaffold_lengths[scaffold_name]
        
        # Get mapping info for this scaffold
        chr_mappings = mapping_stats.get(scaffold_name, {})
        
        # Calculate total aligned bases and percentages for each chromosome
        total_aligned = sum(chr_mappings.values())
        
        result = {
            'scaffold_name': scaffold_name,
            'scaffold_length': scaffold_length,
            'total_aligned_bases': total_aligned,
            'percent_aligned': round(100.0 * total_aligned / scaffold_length, 2) if scaffold_length > 0 else 0.0
        }
        
        # Add per-chromosome mapping percentages
        for chr_name in chromosomes:
            aligned_bases = chr_mappings.get(chr_name, 0)
            pct_aligned = round(100.0 * aligned_bases / scaffold_length, 2) if scaffold_length > 0 else 0.0
            result[f'{chr_name}_aligned_bases'] = aligned_bases
            result[f'{chr_name}_pct'] = pct_aligned
        
        # Find best matching chromosome
        if chr_mappings:
            best_chr = max(chr_mappings.keys(), key=lambda k: chr_mappings[k])
            best_chr_bases = chr_mappings[best_chr]
            best_chr_pct = round(100.0 * best_chr_bases / scaffold_length, 2) if scaffold_length > 0 else 0.0
        else:
            best_chr = "none"
            best_chr_bases = 0
            best_chr_pct = 0.0
        
        result['best_match_chr'] = best_chr
        result['best_match_aligned_bases'] = best_chr_bases
        result['best_match_pct'] = best_chr_pct
        
        results.append(result)
    
    return results, chromosomes


def write_output(results, chromosomes, output_file):
    """
    Write results to TSV file.
    
    Args:
        results: List of dictionaries with mapping statistics
        chromosomes: List of chromosome names
        output_file: Path to output TSV file
    """
    if not results:
        print("Warning: No results to write", file=sys.stderr)
        return
    
    # Define column order
    base_columns = [
        'scaffold_name',
        'scaffold_length',
        'total_aligned_bases',
        'percent_aligned',
        'best_match_chr',
        'best_match_aligned_bases',
        'best_match_pct'
    ]
    
    # Add per-chromosome columns
    chr_columns = []
    for chr_name in chromosomes:
        chr_columns.append(f'{chr_name}_aligned_bases')
        chr_columns.append(f'{chr_name}_pct')
    
    all_columns = base_columns + chr_columns
    
    # Write output
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(all_columns) + '\n')
        
        # Write data
        for result in results:
            values = [str(result.get(col, '')) for col in all_columns]
            f.write('\t'.join(values) + '\n')


def main():
    """Main function."""
    args = parse_arguments()
    
    print(f"Reading scaffold lengths from {args.non_chr_fasta}...", file=sys.stderr)
    scaffold_lengths = get_scaffold_lengths(args.non_chr_fasta)
    print(f"  Found {len(scaffold_lengths)} non-chromosome scaffolds", file=sys.stderr)
    
    print(f"Reading chromosome list from {args.chr_list}...", file=sys.stderr)
    chromosomes = get_chromosome_list(args.chr_list)
    print(f"  Found {len(chromosomes)} chromosomes", file=sys.stderr)
    
    print(f"Parsing BAM file {args.bam}...", file=sys.stderr)
    results, chromosomes = calculate_mapping_stats(
        args.bam,
        scaffold_lengths,
        chromosomes,
        args.min_mapq
    )
    print(f"  Processed {len(results)} scaffolds", file=sys.stderr)
    
    print(f"Writing output to {args.output}...", file=sys.stderr)
    write_output(results, chromosomes, args.output)
    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
