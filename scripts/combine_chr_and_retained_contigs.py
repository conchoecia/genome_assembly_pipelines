#!/usr/bin/env python3
"""
Combine chromosome-scale scaffolds with retained contigs in original assembly order.

This script creates a final curated assembly by:
1. Including all chromosome-scale scaffolds (unchanged)
2. Including retained contigs from non-chromosome scaffolds (after filtering/trimming)
3. Preserving the original order from the input assembly

Usage:
    python combine_chr_and_retained_contigs.py \\
        --input_assembly input.fasta \\
        --chr_list chr_list.txt \\
        --retained_contigs retained.fasta \\
        --contig_map contig_to_scaffold_map.tsv \\
        --output_fasta combined.fasta
"""

import argparse
import sys
from Bio import SeqIO
from collections import defaultdict


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Combine chromosome scaffolds with retained contigs in original order"
    )
    parser.add_argument(
        '--input_assembly',
        required=True,
        help='Original input assembly FASTA file'
    )
    parser.add_argument(
        '--chr_list',
        required=True,
        help='List of chromosome-scale scaffold names (one per line)'
    )
    parser.add_argument(
        '--retained_contigs',
        required=True,
        help='FASTA file with retained contigs (after filtering/trimming)'
    )
    parser.add_argument(
        '--contig_map',
        required=True,
        help='TSV mapping contigs to parent scaffolds'
    )
    parser.add_argument(
        '--output_fasta',
        required=True,
        help='Output combined assembly FASTA file'
    )
    
    return parser.parse_args()


def load_chromosome_list(chr_list_file):
    """Load set of chromosome-scale scaffold names."""
    chr_set = set()
    with open(chr_list_file) as f:
        for line in f:
            line = line.strip()
            if line:
                chr_set.add(line)
    print(f"Loaded {len(chr_set)} chromosome-scale scaffolds", file=sys.stderr)
    return chr_set


def load_contig_to_scaffold_map(contig_map_file):
    """
    Load mapping of contigs to parent scaffolds.
    
    Returns:
        dict: {contig_name: parent_scaffold_name}
    """
    contig_to_scaffold = {}
    with open(contig_map_file) as f:
        # Skip header
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                contig_name = parts[0]
                parent_scaffold = parts[1]
                contig_to_scaffold[contig_name] = parent_scaffold
    
    print(f"Loaded mapping for {len(contig_to_scaffold)} contigs", file=sys.stderr)
    return contig_to_scaffold


def load_retained_contigs(retained_fasta, contig_to_scaffold):
    """
    Load retained contigs and organize by parent scaffold.
    
    Returns:
        dict: {parent_scaffold: [list of SeqRecord objects]}
    """
    scaffold_to_contigs = defaultdict(list)
    retained_count = 0
    
    with open(retained_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            retained_count += 1
            contig_name = record.id
            
            # Handle piece names (e.g., contig_1_piece_1 or contig_1_trimmed)
            # Extract original contig name
            if '_piece_' in contig_name:
                original_contig = contig_name.split('_piece_')[0]
            elif contig_name.endswith('_trimmed'):
                original_contig = contig_name.rsplit('_trimmed', 1)[0]
            else:
                original_contig = contig_name
            
            # Look up parent scaffold
            if original_contig in contig_to_scaffold:
                parent_scaffold = contig_to_scaffold[original_contig]
                scaffold_to_contigs[parent_scaffold].append(record)
            else:
                print(f"Warning: Could not find parent scaffold for {contig_name} (original: {original_contig})", 
                      file=sys.stderr)
    
    print(f"Loaded {retained_count} retained contigs from {len(scaffold_to_contigs)} parent scaffolds", 
          file=sys.stderr)
    return scaffold_to_contigs


def combine_assembly(input_assembly, chr_set, scaffold_to_contigs, output_fasta):
    """
    Combine chromosome scaffolds and retained contigs in original order.
    """
    output_records = []
    chr_count = 0
    retained_contig_count = 0
    skipped_scaffold_count = 0
    
    with open(input_assembly) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            scaffold_name = record.id
            
            if scaffold_name in chr_set:
                # This is a chromosome-scale scaffold - keep it as-is
                output_records.append(record)
                chr_count += 1
            else:
                # This is a non-chromosome scaffold - check for retained contigs
                if scaffold_name in scaffold_to_contigs:
                    # Add all retained contigs from this scaffold
                    contigs = scaffold_to_contigs[scaffold_name]
                    output_records.extend(contigs)
                    retained_contig_count += len(contigs)
                else:
                    # No retained contigs from this scaffold (all were filtered/trimmed away)
                    skipped_scaffold_count += 1
    
    # Write output
    with open(output_fasta, 'w') as out_handle:
        SeqIO.write(output_records, out_handle, "fasta")
    
    print(f"\nSummary:", file=sys.stderr)
    print(f"  Chromosome-scale scaffolds: {chr_count}", file=sys.stderr)
    print(f"  Retained contigs: {retained_contig_count}", file=sys.stderr)
    print(f"  Non-chromosome scaffolds with no retained contigs: {skipped_scaffold_count}", file=sys.stderr)
    print(f"  Total sequences in output: {len(output_records)}", file=sys.stderr)
    print(f"\nWrote combined assembly to {output_fasta}", file=sys.stderr)


def main():
    args = parse_arguments()
    
    # Load data
    chr_set = load_chromosome_list(args.chr_list)
    contig_to_scaffold = load_contig_to_scaffold_map(args.contig_map)
    scaffold_to_contigs = load_retained_contigs(args.retained_contigs, contig_to_scaffold)
    
    # Combine and write output
    combine_assembly(
        args.input_assembly,
        chr_set,
        scaffold_to_contigs,
        args.output_fasta
    )


if __name__ == '__main__':
    main()
