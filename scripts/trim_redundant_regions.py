#!/usr/bin/env python3
"""
Trim redundant regions from contigs based on alignment to chromosome scaffolds.

This script:
1. Identifies redundant blocks (e.g., ≥10 kbp with MAPQ ≥20)
2. Breaks contigs at redundant regions
3. Keeps non-redundant pieces above minimum size
4. Removes highly redundant "swiss cheese" contigs
5. Maintains original contig order in output
"""

import argparse
import pysam
import sys
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Trim redundant regions from contigs"
    )
    parser.add_argument(
        "--input_fasta",
        required=True,
        help="Input FASTA file with contigs"
    )
    parser.add_argument(
        "--bam",
        required=True,
        help="BAM file with contig alignments to chromosomes"
    )
    parser.add_argument(
        "--output_fasta",
        required=True,
        help="Output FASTA file with trimmed contigs"
    )
    parser.add_argument(
        "--report",
        required=True,
        help="Output report file with trimming statistics"
    )
    parser.add_argument(
        "--min_block_size",
        type=int,
        default=10000,
        help="Minimum block size to consider redundant (default: 10000)"
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=20,
        help="Minimum mapping quality for redundant blocks (default: 20)"
    )
    parser.add_argument(
        "--min_keep_piece",
        type=int,
        default=5000,
        help="Minimum size of non-redundant pieces to keep (default: 5000)"
    )
    parser.add_argument(
        "--max_redundant_pct",
        type=float,
        default=80.0,
        help="Remove entire contig if >this %% redundant (default: 80.0)"
    )
    
    return parser.parse_args()


def parse_cigar_blocks(cigar, min_gap=1000):
    """
    Parse CIGAR string and return list of aligned blocks with their query positions.
    Splits on large deletions (≥min_gap) which indicate breaks in alignment.
    
    Returns:
        List of (query_start, query_end, aligned_bases) tuples
    """
    blocks = []
    query_pos = 0
    current_block_start = 0
    current_block_bases = 0
    
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(length)
        
        if op in ['M', '=', 'X']:  # Match/mismatch
            current_block_bases += length
            query_pos += length
        elif op == 'I':  # Insertion (query only)
            current_block_bases += length
            query_pos += length
        elif op == 'D':  # Deletion (reference only)
            if length >= min_gap:  # Large deletion = break in alignment
                if current_block_bases > 0:
                    blocks.append((
                        current_block_start,
                        query_pos,
                        current_block_bases
                    ))
                current_block_start = query_pos
                current_block_bases = 0
        elif op == 'S':  # Soft clip
            query_pos += length
            if current_block_bases > 0:  # End of block
                blocks.append((
                    current_block_start,
                    query_pos - length,  # Don't include soft clip
                    current_block_bases
                ))
            current_block_start = query_pos
            current_block_bases = 0
        elif op == 'H':  # Hard clip (doesn't consume query)
            if current_block_bases > 0:
                blocks.append((
                    current_block_start,
                    query_pos,
                    current_block_bases
                ))
            current_block_bases = 0
    
    # Add final block
    if current_block_bases > 0:
        blocks.append((current_block_start, query_pos, current_block_bases))
    
    return blocks


def identify_redundant_regions(bam_file, min_block_size, min_mapq):
    """
    Identify redundant regions in each contig based on alignment blocks.
    
    Returns:
        Dictionary: {contig_name: [(start, end), ...]} of redundant regions
    """
    redundant_regions = defaultdict(list)
    contig_lengths = {}
    
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    
    for read in bamfile:
        if read.is_unmapped:
            continue
        
        # Skip low quality and supplementary alignments
        if read.mapping_quality < min_mapq:
            continue
        if read.is_supplementary:  # Skip supplementary (flag 2048)
            continue
        
        query_name = read.query_name
        cigar = read.cigarstring
        
        # Track contig lengths
        if query_name not in contig_lengths:
            # Get full query length from alignment
            query_length = read.query_length
            if read.query_alignment_start is not None:
                query_length = read.infer_query_length()
            contig_lengths[query_name] = query_length
        
        # Parse CIGAR to get alignment blocks
        blocks = parse_cigar_blocks(cigar)
        
        # Add blocks that meet size threshold as redundant regions
        for block_start, block_end, block_size in blocks:
            if block_size >= min_block_size:
                redundant_regions[query_name].append((block_start, block_end))
    
    bamfile.close()
    
    # Merge overlapping redundant regions for each contig
    merged_regions = {}
    for contig, regions in redundant_regions.items():
        merged_regions[contig] = merge_intervals(regions)
    
    return merged_regions, contig_lengths


def merge_intervals(intervals):
    """
    Merge overlapping intervals.
    
    Args:
        intervals: List of (start, end) tuples
        
    Returns:
        Sorted list of merged (start, end) tuples
    """
    if not intervals:
        return []
    
    # Sort by start position
    sorted_intervals = sorted(intervals)
    merged = [sorted_intervals[0]]
    
    for start, end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]
        
        if start <= last_end:  # Overlapping or adjacent
            # Merge by extending the last interval
            merged[-1] = (last_start, max(last_end, end))
        else:
            # No overlap, add as new interval
            merged.append((start, end))
    
    return merged


def calculate_non_redundant_pieces(contig_length, redundant_regions, min_piece_size):
    """
    Calculate non-redundant pieces from a contig.
    
    Args:
        contig_length: Total length of contig
        redundant_regions: List of (start, end) redundant intervals
        min_piece_size: Minimum size of pieces to keep
        
    Returns:
        List of (start, end) tuples for non-redundant pieces
    """
    if not redundant_regions:
        # No redundant regions, keep entire contig
        return [(0, contig_length)]
    
    pieces = []
    last_end = 0
    
    for start, end in sorted(redundant_regions):
        # Add non-redundant piece before this redundant region
        if start > last_end:
            piece_size = start - last_end
            if piece_size >= min_piece_size:
                pieces.append((last_end, start))
        
        last_end = max(last_end, end)
    
    # Add final piece after last redundant region
    if last_end < contig_length:
        piece_size = contig_length - last_end
        if piece_size >= min_piece_size:
            pieces.append((last_end, contig_length))
    
    return pieces


def process_contigs(input_fasta, redundant_regions, contig_lengths, 
                   min_piece_size, max_redundant_pct):
    """
    Process contigs and extract non-redundant pieces.
    
    Returns:
        Tuple of (output_records, statistics)
    """
    output_records = []
    stats = {
        'total_contigs': 0,
        'kept_intact': 0,
        'trimmed': 0,
        'removed_swiss_cheese': 0,
        'removed_no_pieces': 0,
        'total_pieces_output': 0,
        'total_bp_input': 0,
        'total_bp_output': 0,
        'total_bp_trimmed': 0
    }
    
    with open(input_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            stats['total_contigs'] += 1
            stats['total_bp_input'] += len(record.seq)
            
            contig_name = record.id
            contig_length = len(record.seq)
            
            # Get redundant regions for this contig
            redun_regions = redundant_regions.get(contig_name, [])
            
            # Calculate redundant percentage
            if redun_regions:
                total_redundant = sum(end - start for start, end in redun_regions)
                redundant_pct = 100.0 * total_redundant / contig_length
            else:
                total_redundant = 0
                redundant_pct = 0.0
            
            # Check if contig is too redundant ("swiss cheese")
            if redundant_pct > max_redundant_pct:
                stats['removed_swiss_cheese'] += 1
                continue
            
            # Calculate non-redundant pieces
            pieces = calculate_non_redundant_pieces(
                contig_length, redun_regions, min_piece_size
            )
            
            if not pieces:
                # No pieces large enough to keep
                stats['removed_no_pieces'] += 1
                continue
            
            # Generate output records
            if len(pieces) == 1 and pieces[0] == (0, contig_length):
                # Keep entire contig intact
                output_records.append(record)
                stats['kept_intact'] += 1
                stats['total_pieces_output'] += 1
                stats['total_bp_output'] += len(record.seq)
            else:
                # Contig was trimmed - output pieces
                stats['trimmed'] += 1
                
                for piece_idx, (start, end) in enumerate(pieces, 1):
                    piece_seq = record.seq[start:end]
                    piece_length = len(piece_seq)
                    
                    # Create piece name
                    if len(pieces) == 1:
                        piece_name = f"{contig_name}_trimmed"
                    else:
                        piece_name = f"{contig_name}_piece_{piece_idx}"
                    
                    piece_record = SeqRecord(
                        piece_seq,
                        id=piece_name,
                        description=f"trimmed from {contig_name} pos {start}-{end}"
                    )
                    
                    output_records.append(piece_record)
                    stats['total_pieces_output'] += 1
                    stats['total_bp_output'] += piece_length
                
                # Track trimmed bases
                stats['total_bp_trimmed'] += (contig_length - sum(end - start for start, end in pieces))
    
    return output_records, stats


def write_report(report_file, stats, args):
    """Write trimming report."""
    with open(report_file, 'w') as f:
        f.write("# Redundant Region Trimming Report\n")
        f.write("#\n")
        f.write("# Parameters:\n")
        f.write(f"#   Min block size: {args.min_block_size:,} bp\n")
        f.write(f"#   Min MAPQ: {args.min_mapq}\n")
        f.write(f"#   Min keep piece: {args.min_keep_piece:,} bp\n")
        f.write(f"#   Max redundant %%: {args.max_redundant_pct:.1f}%%\n")
        f.write("#\n")
        f.write("# Results:\n")
        f.write(f"#   Total input contigs: {stats['total_contigs']}\n")
        f.write(f"#   Kept intact: {stats['kept_intact']}\n")
        f.write(f"#   Trimmed: {stats['trimmed']}\n")
        f.write(f"#   Removed (>80%% redundant): {stats['removed_swiss_cheese']}\n")
        f.write(f"#   Removed (no pieces ≥{args.min_keep_piece} bp): {stats['removed_no_pieces']}\n")
        f.write(f"#   Total output pieces: {stats['total_pieces_output']}\n")
        f.write("#\n")
        f.write(f"#   Total input bp: {stats['total_bp_input']:,}\n")
        f.write(f"#   Total output bp: {stats['total_bp_output']:,}\n")
        f.write(f"#   Total trimmed bp: {stats['total_bp_trimmed']:,}\n")
        f.write(f"#   Reduction: {100.0 * stats['total_bp_trimmed'] / stats['total_bp_input']:.2f}%%\n")


def main():
    """Main execution function."""
    args = parse_arguments()
    
    print(f"Trimming redundant regions from contigs...", file=sys.stderr)
    print(f"  Parameters:", file=sys.stderr)
    print(f"    Min block size: {args.min_block_size:,} bp", file=sys.stderr)
    print(f"    Min MAPQ: {args.min_mapq}", file=sys.stderr)
    print(f"    Min keep piece: {args.min_keep_piece:,} bp", file=sys.stderr)
    print(f"    Max redundant: {args.max_redundant_pct:.1f}%", file=sys.stderr)
    print(file=sys.stderr)
    
    # Step 1: Identify redundant regions
    print("Step 1: Identifying redundant regions from BAM...", file=sys.stderr)
    redundant_regions, contig_lengths = identify_redundant_regions(
        args.bam,
        args.min_block_size,
        args.min_mapq
    )
    print(f"  Found redundant regions in {len(redundant_regions)} contigs", 
          file=sys.stderr)
    print(file=sys.stderr)
    
    # Step 2: Process contigs and extract non-redundant pieces
    print("Step 2: Processing contigs and extracting non-redundant pieces...", 
          file=sys.stderr)
    output_records, stats = process_contigs(
        args.input_fasta,
        redundant_regions,
        contig_lengths,
        args.min_keep_piece,
        args.max_redundant_pct
    )
    print(file=sys.stderr)
    
    # Step 3: Write output
    print("Step 3: Writing output...", file=sys.stderr)
    with open(args.output_fasta, 'w') as out_handle:
        SeqIO.write(output_records, out_handle, "fasta")
    
    write_report(args.report, stats, args)
    
    # Print summary
    print(f"\nTrimming complete:", file=sys.stderr)
    print(f"  Input contigs: {stats['total_contigs']}", file=sys.stderr)
    print(f"  Kept intact: {stats['kept_intact']}", file=sys.stderr)
    print(f"  Trimmed: {stats['trimmed']}", file=sys.stderr)
    print(f"  Removed (swiss cheese): {stats['removed_swiss_cheese']}", file=sys.stderr)
    print(f"  Removed (too small): {stats['removed_no_pieces']}", file=sys.stderr)
    print(f"  Total output pieces: {stats['total_pieces_output']}", file=sys.stderr)
    print(f"  Input bp: {stats['total_bp_input']:,}", file=sys.stderr)
    print(f"  Output bp: {stats['total_bp_output']:,}", file=sys.stderr)
    print(f"  Trimmed: {stats['total_bp_trimmed']:,} bp ({100.0 * stats['total_bp_trimmed'] / stats['total_bp_input']:.2f}%)", 
          file=sys.stderr)
    print(f"\n  Output written to: {args.output_fasta}", file=sys.stderr)
    print(f"  Report written to: {args.report}", file=sys.stderr)


if __name__ == '__main__':
    main()
