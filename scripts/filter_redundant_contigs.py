#!/usr/bin/env python3
"""
Filter redundant contigs based on alignment coverage and mapping quality.

Removes contigs that have high-quality alignments to chromosome-scale scaffolds,
indicating they are redundant with the chromosome assembly.
"""

import argparse
import pandas as pd
from Bio import SeqIO
import sys


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter redundant contigs from assembly"
    )
    parser.add_argument(
        "--input_fasta",
        required=True,
        help="Input FASTA file with contigs"
    )
    parser.add_argument(
        "--summary_csv",
        required=True,
        help="CSV file with contig classification summary"
    )
    parser.add_argument(
        "--output_fasta",
        required=True,
        help="Output FASTA file with filtered contigs"
    )
    parser.add_argument(
        "--removed_list",
        required=True,
        help="Output file listing removed contigs"
    )
    parser.add_argument(
        "--min_coverage",
        type=float,
        default=95.0,
        help="Minimum coverage percentage for filtering (default: 95.0)"
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=30,
        help="Minimum mapping quality for filtering (default: 30)"
    )
    parser.add_argument(
        "--min_lr_coverage",
        type=float,
        default=None,
        help="Optional: Maximum long-read coverage for low-coverage filtering"
    )
    
    return parser.parse_args()


def identify_contigs_to_remove(summary_csv, min_coverage, min_mapq, min_lr_coverage=None):
    """
    Identify contigs to remove based on filtering criteria.
    
    Args:
        summary_csv: Path to classification summary CSV
        min_coverage: Minimum alignment coverage percentage
        min_mapq: Minimum mapping quality
        min_lr_coverage: Optional maximum long-read coverage threshold
        
    Returns:
        Tuple of (set of contig names to remove, DataFrame with statistics)
    """
    # Read summary
    df = pd.read_csv(summary_csv)
    
    # Ensure required columns exist
    required_cols = ['contig_name', 'best_match_pct', 'best_match_mapq']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' not found in summary CSV")
    
    # Apply filters
    remove_mask = (df['best_match_pct'] >= min_coverage) & (df['best_match_mapq'] >= min_mapq)
    
    # Optional: also filter low coverage contigs
    if min_lr_coverage is not None and 'mean_coverage' in df.columns:
        low_coverage_mask = df['mean_coverage'] < min_lr_coverage
        remove_mask = remove_mask | low_coverage_mask
    
    # Get contigs to remove
    contigs_to_remove = set(df[remove_mask]['contig_name'].values)
    
    # Get statistics
    removed_df = df[remove_mask].copy()
    kept_df = df[~remove_mask].copy()
    
    return contigs_to_remove, removed_df, kept_df


def filter_fasta(input_fasta, output_fasta, contigs_to_remove):
    """
    Filter FASTA file, removing specified contigs.
    
    Args:
        input_fasta: Input FASTA file path
        output_fasta: Output FASTA file path
        contigs_to_remove: Set of contig names to exclude
        
    Returns:
        Tuple of (num_kept, num_removed, total_bp_kept, total_bp_removed)
    """
    num_kept = 0
    num_removed = 0
    total_bp_kept = 0
    total_bp_removed = 0
    
    kept_records = []
    
    with open(input_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in contigs_to_remove:
                num_removed += 1
                total_bp_removed += len(record.seq)
            else:
                kept_records.append(record)
                num_kept += 1
                total_bp_kept += len(record.seq)
    
    # Write filtered FASTA
    with open(output_fasta, 'w') as out_handle:
        SeqIO.write(kept_records, out_handle, "fasta")
    
    return num_kept, num_removed, total_bp_kept, total_bp_removed


def write_removed_list(removed_df, removed_list_path, filter_stats):
    """
    Write list of removed contigs with statistics.
    
    Args:
        removed_df: DataFrame with removed contigs
        removed_list_path: Output file path
        filter_stats: Tuple of (num_kept, num_removed, total_bp_kept, total_bp_removed)
    """
    num_kept, num_removed, total_bp_kept, total_bp_removed = filter_stats
    
    with open(removed_list_path, 'w') as f:
        # Write summary
        f.write("# Redundant Contig Filtering Results\n")
        f.write(f"# Contigs kept: {num_kept}\n")
        f.write(f"# Contigs removed: {num_removed}\n")
        f.write(f"# Total bp kept: {total_bp_kept:,}\n")
        f.write(f"# Total bp removed: {total_bp_removed:,}\n")
        f.write("#\n")
        
        # Write column headers
        cols_to_write = ['contig_name', 'contig_length', 'best_match_chr', 
                        'best_match_pct', 'best_match_mapq']
        if 'mean_coverage' in removed_df.columns:
            cols_to_write.append('mean_coverage')
        
        f.write('\t'.join(cols_to_write) + '\n')
        
        # Write removed contigs
        for _, row in removed_df.iterrows():
            values = [str(row[col]) for col in cols_to_write]
            f.write('\t'.join(values) + '\n')


def main():
    """Main execution function."""
    args = parse_arguments()
    
    print(f"Filtering redundant contigs...", file=sys.stderr)
    print(f"  Criteria: coverage >= {args.min_coverage}%, MAPQ >= {args.min_mapq}", 
          file=sys.stderr)
    
    # Identify contigs to remove
    contigs_to_remove, removed_df, kept_df = identify_contigs_to_remove(
        args.summary_csv,
        args.min_coverage,
        args.min_mapq,
        args.min_lr_coverage
    )
    
    print(f"  Found {len(contigs_to_remove)} contigs to remove", file=sys.stderr)
    
    # Filter FASTA
    filter_stats = filter_fasta(
        args.input_fasta,
        args.output_fasta,
        contigs_to_remove
    )
    
    num_kept, num_removed, total_bp_kept, total_bp_removed = filter_stats
    
    # Write removed list
    write_removed_list(removed_df, args.removed_list, filter_stats)
    
    # Print summary
    print(f"\nFiltering complete:", file=sys.stderr)
    print(f"  Contigs kept: {num_kept} ({total_bp_kept:,} bp)", file=sys.stderr)
    print(f"  Contigs removed: {num_removed} ({total_bp_removed:,} bp)", file=sys.stderr)
    print(f"  Filtered assembly written to: {args.output_fasta}", file=sys.stderr)
    print(f"  Removed contigs list: {args.removed_list}", file=sys.stderr)


if __name__ == '__main__':
    main()
