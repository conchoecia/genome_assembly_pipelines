#!/usr/bin/env python3
"""
This program prints out the position of any given kmers. Used with bedtools merge it is 
useful to get ranges and sizes of telomere or other repeat stretches.

Use like this:
cat genome.fasta | python kmer_positions.py TTAGGG TTAAGG --canonical | bedtools merge | awk '{if (($3-$2)/6 > 1) {printf("%s\t%d\t%d\t%d\n", $1, $2, $3, ($3-$2)/6 )}}'

Flags:
  --canonical : match both the input kmer(s) and their reverse complements
"""

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("kmers", nargs='+', help="Kmers to search for (must be same length)")
    parser.add_argument("--canonical", action="store_true",
                        help="Also match reverse complements of kmers")
    return parser.parse_args()

args = parse_args()

# Validate and prepare kmers
input_kmers = [k.upper() for k in args.kmers]
kmer_sizes = set(len(k) for k in input_kmers)
if len(kmer_sizes) != 1:
    raise ValueError(f"All kmers must be the same length. Got: {sorted(kmer_sizes)}")

k = list(kmer_sizes)[0]
tkmers = set(input_kmers)
if args.canonical:
    rc_kmers = set(reverse_complement(k) for k in input_kmers)
    tkmers.update(rc_kmers)

# Search genome
for record in SeqIO.parse(sys.stdin, "fasta"):
    seq = str(record.seq).upper()
    for i in range(len(seq) - k + 1):
        if seq[i:i+k] in tkmers:
            print("\t".join([record.id, str(i), str(i + k)]))