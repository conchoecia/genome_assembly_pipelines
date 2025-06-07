#!/usr/bin/env python

"""
Title:  NCBI_fcs_gx_filter.py
Author: Darrin T. Schultz https://github.com/conchoecia
Date:   May 2025

This script takes a FCS-gx report file from the NCBI CLI tool
 and splits the genome where the contaminating sequences were.

This tool is designed to identify, for example, prokaryotic sequence
 contamination in a eukaryotic genome assembly.

To be used with NCBI FCS: https://github.com/ncbi/fcs
Specifically, with the gx tool: https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart
"""

import argparse
import os
import pandas
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# parse the args. We need a fasta file and the contamination file.
def parse_args():
    """
    Parse the arguments

    We need the fasta file and the contamination file.
    Print help if no arguments are given
    """
    parser = argparse.ArgumentParser(description="Split or mask a genome at the contamination sites")
    parser.add_argument("-f", "--fasta_input", help="The fasta file to process")
    parser.add_argument("-r", "--report", required = True,
                        help="The report file from NCBI FCS-gx")
    parser.add_argument("-k", "--keep_fasta",   help="This is the fasta file of the KEPT (non-contam) sequences.", required=True)
    parser.add_argument("-c", "--contam_fasta", help="This is the fasta file of the REJECTED (contam) sequences.", required=True)
    parser.add_argument("-n", "--mask-with-n", action="store_true", help="Mask contaminant regions with Ns instead of splitting")
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.fasta_input):
        raise IOError("The fasta file does not exist")
    if not os.path.isfile(args.report):
        raise IOError("The report file does not exist")

    # ensure that the keep_fasta, contam_fasta, and the fasta are not the same file
    if len(list(set([args.fasta_input, args.keep_fasta, args.contam_fasta]))) < 3:
        raise IOError("The fasta_input, keep_fasta, and contam_fasta files must be different")

    return args

def parse_fcs_gx_contamination_txt(filepath):
    """
    Parse FCS-GX-style contamination report to extract regions to exclude.

    Expected format:
    - Header lines starting with "#" or metadata JSON on the first line
    - Tab-delimited data lines with columns:
      seq_id, start_pos, end_pos, seq_len, action, div, agg_cont_cov, top_tax_name

    Returns:
        A dict mapping sequence IDs to sets of (start, stop) tuples (0-based, exclusive end).
    """
    trim_regions = {}

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("["):
                continue  # skip headers and metadata

            parts = line.split("\t")
            if len(parts) < 3 or parts[4] not in ["EXCLUDE", "REVIEW"]:
                if parts[4] not in ["EXCLUDE", "REVIEW"]:
                    print(parts)
                    raise IOError("We only know how to handle EXCLUDE and REVIEW lines.")
                continue  # skip malformed or non-EXCLUDE lines

            seq_id = parts[0]
            start = int(parts[1]) - 1  # convert to 0-based
            stop = int(parts[2])       # exclusive end

            if seq_id not in trim_regions:
                trim_regions[seq_id] = set()
            trim_regions[seq_id].add((start, stop))

    return trim_regions

def parse_fasta(fasta_input, keep_fasta, contam_fasta,
                trim_regions, mask_with_n=False):
    """
    Process a FASTA file by removing or masking contaminant regions.

    Args:
        fasta_input (str): Path to input FASTA file.
        keep_fasta (str): Path to output FASTA for retained sequences.
        contam_fasta (str): Path to output FASTA for removed contaminant sequences.
        trim_regions (dict): Dict mapping sequence IDs to sets of (start, stop) tuples.
        mask_with_n (bool): If True, contaminant regions are masked with Ns instead of being removed.
    """
    handle = open(fasta_input, "r")
    keep_out = open(keep_fasta, "w")
    contam_out = open(contam_fasta, "w")

    for record in SeqIO.parse(handle, "fasta"):
        scaf_counter = 1
        removed_counter = 1

        if record.id not in trim_regions:
            SeqIO.write(record, keep_out, "fasta")
            continue

        if mask_with_n:
            new_seq = list(str(record.seq))
            for start, stop in trim_regions[record.id]:
                for i in range(start, stop):  # stop is exclusive
                    if i < len(new_seq):
                        new_seq[i] = 'N'
            masked_record = SeqRecord(Seq("".join(new_seq)), id=record.id, description=record.description)
            SeqIO.write(masked_record, keep_out, "fasta")
        else:
            ranges_to_remove = trim_regions[record.id]
            keep_ranges = invert_ranges(ranges_to_remove, len(record.seq))

            print(f"Keeping these ranges in {record.id}:", file=sys.stderr)
            print(keep_ranges, file=sys.stderr)

            # Write retained fragments
            for keep_seq in extract_ranges(str(record.seq), keep_ranges):
                if len(keep_seq) < 10:
                    continue
                seqpiece = f"{record.id}::piece{scaf_counter}"
                sr = SeqRecord(Seq(keep_seq), seqpiece, '', '')
                SeqIO.write(sr, keep_out, "fasta")
                scaf_counter += 1

            # Write removed fragments
            for trim_seq in extract_ranges(str(record.seq), sorted(ranges_to_remove)):
                if len(trim_seq) < 10:
                    continue
                seqpiece = f"{record.id}::removed{removed_counter}"
                sr = SeqRecord(Seq(trim_seq), seqpiece, '', '')
                SeqIO.write(sr, contam_out, "fasta")
                removed_counter += 1

    handle.close()
    keep_out.close()
    contam_out.close()

def invert_ranges(ranges_to_remove, string_length):
    """
    # Example usage:
    string_length = 20
    ranges_to_remove = [(2, 5), (8, 11), (14, 17)]
    ranges_to_keep = invert_ranges(ranges_to_remove, string_length)
    print(ranges_to_keep)  # Output: [(0, 1), (6, 7), (12, 13), (18, 19)]
    """
    ranges_to_keep = []
    start = 0

    for range_start, range_end in sorted(ranges_to_remove):
        if range_start > string_length or range_end < 0:
            continue

        if range_start > start:
            # we don't used range_start - 1 because of 0-based indexing
            #  and the face that python slicing goes up until this index
            ranges_to_keep.append((start, range_start))

        start = max(start, range_end + 1)

    if start < string_length:
        # Same as above. We don't do string_length-1 because of python slicing.
        ranges_to_keep.append((start, string_length))

    return ranges_to_keep

def extract_ranges(input_string, ranges_to_keep):
    for range_start, range_end in sorted(ranges_to_keep):
        if range_start >= len(input_string) or range_end < 0:
            raise ValueError("Invalid range: ({}, {}). The range is outside the bounds of the input string.".format(range_start, range_end))

        yield input_string[range_start:range_end]

def print_trim_regions(trim_regions):
    """
    Print the trim regions for each scaffold to stderr in a readable format.
    """
    for chrom in sorted(trim_regions):
        print(f"{chrom}:", file=sys.stderr)
        for start, stop in sorted(trim_regions[chrom]):
            print(f"  Trim region: {start}..{stop}", file=sys.stderr)

def main():
    args = parse_args()

    trim_regions = parse_fcs_gx_contamination_txt(args.report)
    print_trim_regions(trim_regions)  # optional
    parse_fasta(args.fasta_input, args.keep_fasta, args.contam_fasta,
                trim_regions, mask_with_n=args.mask_with_n)

if __name__ == "__main__":
    main()