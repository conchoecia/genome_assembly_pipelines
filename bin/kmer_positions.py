#!/usr/bin/env python3
"""
This program prints out the position of any given kmers. Used with bedtools merge it is 
useful to get ranges and sizes of telomere or other repeat stretches.

first and later positional args are kmers. They all must be the same size.

use like this:
cat genome.fasta | python kmer_positions.py TTAGGG TTAAGG | bedtools merge | awk '{if (($3-$2)/6 > 1) {printf("%s\t%d\t%d\t%d\n", $1, $2, $3, ($3-$2)/6 )}}'

^ prints out a bedgraph file (incomplete, no 0 stretches), showing how many kmers are in stretches.
"""

from Bio import SeqIO
import sys

tkmers = list(set([x.upper() for x in sys.argv[1::]]))
kmer_sizes = list(set([len(x) for x in tkmers]))
if len(kmer_sizes) > 1:
    raise IOError("All kmers must have the same length. We saw k = {}".format(kmer_sizes))
k = kmer_sizes[0]

for record in SeqIO.parse(sys.stdin, "fasta"):
    positions = []
    thisseq = str(record.seq)
    #while not done:
    for i in range(len(thisseq)):
        if thisseq[i:i+k].upper() in tkmers:
            print("\t".join([record.id, str(i), str(i+k)]))
