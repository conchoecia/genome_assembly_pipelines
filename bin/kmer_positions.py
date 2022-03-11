#!/usr/bin/env python3
"""
This program prints out the position of any given kmers. Used with bedtools merge it is 
useful to get ranges and sizes of telomere or other repeat stretches.

first positional arg is k size
second and later positional args are kmers

use like this:
cat genome.fasta | python telomere_peaks.py 6 TTAGGG CCCTAA | bedtools merge | awk '{if (($3-$2)/6 > 1) {printf("%s\t%d\t%d\t%d\n", $1, $2, $3, ($3-$2)/6 )}}'

^ prints out a bedgraph file (incomplete, no 0 stretches), showing how many kmers are in stretches.
"""

from Bio import SeqIO
import sys

def reverse_complement(seq):
    this_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(this_complement.get(base, base) for base in reversed(seq))

assembly = "/bigdata/user/darrin/pyrosome/assembly_versions/pyro_hapA_v0.1.fasta"

k = int(sys.argv[1])
kmers = [x.upper() for x in sys.argv[2::]]
revcomp = [reverse_complement(x) for x in kmers]
tkmers = list(set(kmers + revcomp))

for record in SeqIO.parse(sys.stdin, "fasta"):
    positions = []
    thisseq = str(record.seq)
    #while not done:
    for i in range(len(thisseq)):
        if thisseq[i:i+k].upper() in tkmers:
            print("\t".join([record.id, str(i), str(i+k)]))
