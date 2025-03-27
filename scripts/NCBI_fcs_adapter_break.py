#!/usr/bin/env python

"""
Title:  NCBI_fcs_adapter_break.py
Author: Darrin T. Schultz https://github.com/conchoecia
Date:   July 2023

This script takes an adapter contamination file from NCBI FCS-adaptor and splits a genome
where the vector contamination was.

To be used with NCBI FCS-adaptor: https://github.com/ncbi/fcs
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
    parser.add_argument("-f", "--fasta", help="The fasta file to process")
    parser.add_argument("-c", "--contamination", help="The contamination file from NCBI FCS-adaptor")
    parser.add_argument("--ncbi", action="store_true",
                    help="Use NCBI RemainingContamination.txt format instead of FCS-adaptor contamination file")
    parser.add_argument("-n", "--mask-with-n", action="store_true", help="Mask contaminant regions with Ns instead of splitting")
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.fasta):
        raise IOError("The fasta file does not exist")
    if not os.path.isfile(args.contamination):
        raise IOError("The contamination file does not exist")

    return args

def parse_ncbi_remaining_contamination_txt(contam_filepath):
    """
    Parse NCBI RemainingContamination.txt-style file to extract trim regions.

    Example lines:
    Scucuzx_259;HRSCAF=7880	91568	89088..89190	vector/etc-not_cleaned
    """
    trim_regions = {}
    in_trim_section = False

    with open(contam_filepath, "r") as f:
        for line in f:
            line = line.strip()

            # Look for start of the "Trim:" section
            if line.startswith("Trim:"):
                in_trim_section = True
                continue

            # Skip until we're inside the Trim section
            if not in_trim_section:
                continue

            # Skip header line in Trim section
            if line.startswith("Sequence name") or not line:
                continue

            # Parse the scaffold ID and the trim region
            parts = line.split("\t")
            if len(parts) < 3:
                continue

            scaffold_id = parts[0]  # keep full ID
            span_string = parts[2]

            for span in span_string.split(","):
                if ".." in span:
                    start = int(span.split("..")[0]) - 1
                    stop  = int(span.split("..")[1]) - 1
                    if scaffold_id not in trim_regions:
                        trim_regions[scaffold_id] = set()
                    trim_regions[scaffold_id].add((start, stop))

    return trim_regions

def parse_contamination_file(contam_filepath):
    """
    We need to parse the contamination file from NCBI FCS-adaptor

    Return a data structure of the things to keep (or the things to get rid of)
    """
    # open the contamination
    df = pandas.read_csv(contam_filepath, sep="\t")
    df.columns = [x.replace("#","") for x in df.columns]
    # read the contamination file
    trim_regions = {}

    # iterate through the rows of the dataframe
    for index, row in df.iterrows():
        # if the action is to trim, then record this section
        if row["action"] == "ACTION_TRIM":
            chrom = row["accession"]
            ranges = row["range"].split(",")
            for startstop in ranges:
                # subtract 1 from the start because the strings of SeqIO are 0-based
                start = int(startstop.split("..")[0]) - 1
                stop  = int(startstop.split("..")[1]) - 1
                if chrom not in trim_regions:
                    trim_regions[chrom] = set()
                trim_regions[chrom].add((start, stop))
        # OTHER ACTIONS
        else:
            raise ValueError("We only know how to handle ACTION_TRIM.")
    return trim_regions

def parse_fasta(fasta_filepath, trim_regions, mask_with_n=False):
    """
    Read in the fasta file and print out only the parts that
     haven't been marked to trim.

    We don't trust the joins that were previously made by the "trim" sequences,
     therefore we don't hardmask.

    Print to std.out
    """
    outhandle = sys.stdout

    # now read in the next sequences and trim out the bad parts
    with open(fasta_filepath, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            scaf_counter = 1
            if record.id not in trim_regions:
                SeqIO.write(record, outhandle, "fasta")
            else:
                if mask_with_n:
                    new_seq = list(str(record.seq))
                    for start, stop in trim_regions[record.id]:
                        for i in range(start, stop + 1):
                            if i < len(new_seq):
                                new_seq[i] = 'N'
                    masked_record = SeqRecord(Seq("".join(new_seq)), id=record.id, description=record.description)
                    SeqIO.write(masked_record, outhandle, "fasta")
                else:
                    ranges_to_remove = trim_regions[record.id]
                    keep_ranges = invert_ranges(ranges_to_remove, len(record.seq))
                    print("Keeping these ranges in {}:".format(record.id), file=sys.stderr)
                    print(keep_ranges, file=sys.stderr)
                    for keep_seq in extract_ranges(str(record.seq), keep_ranges):
                        seqpiece = record.id + "::" + "piece" + str(scaf_counter)
                        sr = SeqRecord(Seq(keep_seq), seqpiece, '', '')
                        SeqIO.write(sr, outhandle, "fasta")
                        scaf_counter += 1

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

        yield input_string[range_start:range_end + 1]

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

    if args.ncbi:
        trim_regions = parse_ncbi_remaining_contamination_txt(args.contamination)
    else:
        trim_regions = parse_contamination_file(args.contamination)

    print_trim_regions(trim_regions)  # optional
    parse_fasta(args.fasta, trim_regions, mask_with_n=args.mask_with_n)


if __name__ == "__main__":
    main()