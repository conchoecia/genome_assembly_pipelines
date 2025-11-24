#!/usr/bin/env python3
"""
Bin Hi-C pairs into a contact matrix.

Reads filtered pairs.gz files and bins contacts into a 2D matrix with configurable
bin size. Stores only the upper triangle for efficiency (matching pairs.gz format).

Usage:
    python bin_pairs.py --input pairs1.gz [pairs2.gz ...] --output matrix.npz --bin-size 2500000
"""

import argparse
import gzip
import sys
from collections import defaultdict
import numpy as np
import scipy.sparse

def parse_header(pairs_file):
    """
    Parse chromosome information from pairs.gz header.
    
    Returns:
        list: ordered list of (chromosome name, length) tuples
    """
    chroms = []
    with gzip.open(pairs_file, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                break
            if line.startswith('#chromosome:'):
                # Format: #chromosome: NAME LENGTH
                parts = line.strip().split()
                if len(parts) >= 3:
                    chrom_name = parts[1]
                    chrom_length = int(parts[2])
                    chroms.append((chrom_name, chrom_length))
    return chroms

def create_diploid_chromosome_map(chroms, bin_size):
    """
    Create mapping from (chrom, phase) to bin indices.
    
    Args:
        chroms: ordered list of (chromosome name, length) tuples
        bin_size: size of bins in bp
    
    Returns:
        tuple: (chrom_map, chrom_list, total_bins)
            - chrom_map: dict of (chrom, phase) -> (start_bin, num_bins)
            - chrom_list: list of (chrom, phase) tuples in order
            - total_bins: total number of bins across all chromosomes
    """
    chrom_map = {}
    chrom_list = []
    current_bin = 0
    
    print("\nChromosome bin allocation:", file=sys.stderr)
    
    # Create diploid chromosomes: chr1.0, chr1.1, chr2.0, chr2.1, etc.
    # Iterate over list to preserve order from pairs.gz header
    for chrom_name, chrom_length in chroms:
        num_bins = (chrom_length + bin_size - 1) // bin_size  # Ceiling division
        
        for phase in [0, 1]:
            chrom_map[(chrom_name, phase)] = (current_bin, num_bins)
            chrom_list.append((chrom_name, phase))
            print(f"  {chrom_name}.{phase}: bins {current_bin}-{current_bin + num_bins - 1} ({num_bins} bins, {chrom_length:,} bp)", 
                  file=sys.stderr)
            current_bin += num_bins
    
    print(f"\nTotal bins allocated: {current_bin:,}\n", file=sys.stderr)
    
    return chrom_map, chrom_list, current_bin

def bin_pairs(pairs_files, chrom_map, bin_size):
    """
    Bin pairs into contact matrix (upper triangle only).
    
    Args:
        pairs_files: list of pairs.gz file paths
        chrom_map: mapping from (chrom, phase) -> (start_bin, num_bins)
        bin_size: size of bins in bp
    
    Returns:
        dict: (bin1, bin2) -> count, where bin1 <= bin2
    """
    contacts = defaultdict(int)
    total_pairs = 0
    skipped_pairs = 0
    
    for pairs_file in pairs_files:
        print(f"Processing {pairs_file}...", file=sys.stderr)
        
        with gzip.open(pairs_file, 'rt') as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                # Parse pairs line
                # Format: readID chr1 pos1 chr2 pos2 strand1 strand2 phase1 phase2
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                chr1, pos1_str, chr2, pos2_str = fields[1], fields[2], fields[3], fields[4]
                phase1_str, phase2_str = fields[7], fields[8]
                
                # Skip unphaseable contacts
                if phase1_str == '.' or phase2_str == '.':
                    skipped_pairs += 1
                    continue
                
                try:
                    pos1 = int(pos1_str)
                    pos2 = int(pos2_str)
                    phase1 = int(phase1_str)
                    phase2 = int(phase2_str)
                except ValueError:
                    skipped_pairs += 1
                    continue
                
                # Get bin indices for each contact end
                if (chr1, phase1) not in chrom_map or (chr2, phase2) not in chrom_map:
                    skipped_pairs += 1
                    continue
                
                start_bin1, _ = chrom_map[(chr1, phase1)]
                start_bin2, _ = chrom_map[(chr2, phase2)]
                
                bin1 = start_bin1 + (pos1 // bin_size)
                bin2 = start_bin2 + (pos2 // bin_size)
                
                # Store upper triangle only (bin1 <= bin2)
                if bin1 <= bin2:
                    contacts[(bin1, bin2)] += 1
                else:
                    contacts[(bin2, bin1)] += 1
                
                total_pairs += 1
                
                if total_pairs % 1000000 == 0:
                    print(f"  Processed {total_pairs:,} pairs, {len(contacts):,} unique bins", 
                          file=sys.stderr)
    
    print(f"Total pairs processed: {total_pairs:,}", file=sys.stderr)
    print(f"Skipped pairs: {skipped_pairs:,}", file=sys.stderr)
    print(f"Non-zero bins: {len(contacts):,}", file=sys.stderr)
    
    return contacts

def save_matrix(contacts, chrom_list, chrom_map, total_bins, bin_size, output_file):
    """
    Save contact matrix as sparse npz file.
    
    Args:
        contacts: dict of (bin1, bin2) -> count
        chrom_list: list of (chrom, phase) tuples
        chrom_map: dict of (chrom, phase) -> (start_bin, num_bins)
        total_bins: total number of bins
        bin_size: size of bins in bp
        output_file: output .npz file path
    """
    print(f"Creating sparse matrix ({total_bins} x {total_bins})...", file=sys.stderr)
    
    # Convert dict to COO format (coordinate format)
    rows = []
    cols = []
    data = []
    
    for (bin1, bin2), count in contacts.items():
        rows.append(bin1)
        cols.append(bin2)
        data.append(count)
    
    # Create sparse matrix
    matrix = scipy.sparse.coo_matrix(
        (data, (rows, cols)),
        shape=(total_bins, total_bins),
        dtype=np.int32
    )
    
    # Convert to CSR for efficient storage
    matrix = matrix.tocsr()
    
    print(f"Matrix density: {matrix.nnz / (total_bins ** 2) * 100:.6f}%", file=sys.stderr)
    print(f"Saving to {output_file}...", file=sys.stderr)
    
    # Save matrix and metadata
    chrom_names = [f"{chrom}.{phase}" for chrom, phase in chrom_list]
    
    scipy.sparse.save_npz(output_file, matrix)
    
    # Save metadata separately with chromosome bin ranges
    metadata_file = output_file.replace('.npz', '.metadata.npz')
    
    # Extract bin ranges for each chromosome from chrom_map
    chrom_start_bins = []
    chrom_num_bins = []
    for chrom, phase in chrom_list:
        start_bin, num_bins = chrom_map[(chrom, phase)]
        chrom_start_bins.append(start_bin)
        chrom_num_bins.append(num_bins)
    
    np.savez(metadata_file,
             chrom_names=chrom_names,
             chrom_list=chrom_list,
             chrom_start_bins=np.array(chrom_start_bins),
             chrom_num_bins=np.array(chrom_num_bins),
             bin_size=bin_size,
             total_bins=total_bins)
    
    print(f"Saved matrix to {output_file}", file=sys.stderr)
    print(f"Saved metadata to {metadata_file}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Bin Hi-C pairs into contact matrix')
    parser.add_argument('--input', '-i', nargs='+', required=True,
                        help='Input pairs.gz file(s)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output .npz file')
    parser.add_argument('--bin-size', '-b', type=int, default=2500000,
                        help='Bin size in bp (default: 2.5 Mb)')
    
    args = parser.parse_args()
    
    # Parse chromosome information from first file
    print("Parsing chromosome information...", file=sys.stderr)
    chroms = parse_header(args.input[0])
    print(f"Found {len(chroms)} chromosomes", file=sys.stderr)
    
    # Create diploid chromosome map
    chrom_map, chrom_list, total_bins = create_diploid_chromosome_map(chroms, args.bin_size)
    print(f"Total bins: {total_bins:,}", file=sys.stderr)
    print(f"Bin size: {args.bin_size:,} bp ({args.bin_size / 1e6:.2f} Mb)", file=sys.stderr)
    
    # Bin pairs
    contacts = bin_pairs(args.input, chrom_map, args.bin_size)
    
    # Save matrix
    save_matrix(contacts, chrom_list, chrom_map, total_bins, args.bin_size, args.output)
    
    print("Done!", file=sys.stderr)

if __name__ == '__main__':
    main()
