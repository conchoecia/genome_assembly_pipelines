#!/usr/bin/env python

"""
This script takes a list of read names, and a list of fastq files.
The script opens the fastq files, and writes those reads to stdout.
"""

#use biopython
import argparse
from Bio import SeqIO
import os
import sys

def parse_args():
    """
    We need to get:
      - a file with the list of read names to keep
      - a list of fastq files. Make this comma-separated
    """
    parser = argparse.ArgumentParser(description="Get reads from a list of fastq files")
    parser.add_argument("-k", "--keepread_file", help="File with the list of reads to keep")
    parser.add_argument("-f", "--fastq", help="Comma-separated list of fastq files")

    args = parser.parse_args()
    # if keepread_file is empty or fastq is empty, print help
    if args.keepread_file is None or args.fastq is None:
        parser.print_help()
        sys.exit(1)
    # for the keep_these_reads, check that the file still exists
    if not os.path.exists(args.keepread_file):
        print("Error: file {} does not exist".format(args.keepread_file))
        sys.exit(1)

    # split the list of fastq files using commas
    fastqs = args.fastq.split(",")
    # check that each one exists
    for fastq in fastqs:
        if not os.path.exists(fastq):
            print("Error: file {} does not exist".format(fastq))
            sys.exit(1)
    args.fastq = fastqs

    return args

def reads_to_keep(readlist_file) -> list:
    """
    Opens the file with the list of reads to keep,
      and returns a list of the read names.
    """
    reads = []
    with open(readlist_file, "r") as f:
        reads = f.readlines()
    return reads

def main():
    # parse the args
    args = parse_args()
    reads = reads_to_keep(args.keepread_file)

    #iterate through the fastq files
    for fastq in args.fastq:
        for record in SeqIO.parse(fastq, "fastq"):
            if record.id in reads:
                print(record.format("fastq").strip())

if __name__ == "__main__":
    main()