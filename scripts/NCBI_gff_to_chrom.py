#!/usr/bin/env python3
"""
Takes an NCBI GFF, fasta, and protein file. Outputs cleaned-up files
for use with ODP.
"""

import os
import sys
from Bio import SeqIO
import pandas as pd

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-a", "--assembly",
                    type=argparse.FileType('r'),
                    help="the path to the genome assembly fasta file.")
parser.add_argument("-g", "--gff",
                    type=argparse.FileType('r'),
                    help="the path to the gff file")
parser.add_argument("-p", "--proteins",
                    type=argparse.FileType('r'),
                    help="the path to the proteins")
parser.add_argument("-r", "--rename",
                    type=argparse.FileType('r'),
                    required = False,
                    help="A path to a two-column tsv file. First col orig fasta header, second col new name")
parser.add_argument("-o", "--output",
                    type=str,
                    required=True,
                    help="the prefix for the output")

# parse the arguments
args = parser.parse_args()

# original script starts here
chrom2loc = os.path.abspath("{}.chrom".format(args.output))
assembly  = os.path.abspath("{}.fasta".format(args.output))
proteins  = os.path.abspath("{}.pep".format(args.output))

#check that these three files don't already exist
files_to_check = [chrom2loc, assembly, proteins]
for f in files_to_check:
    if os.path.exists(f):
        print("ERROR: {} already exists! Rename the output -o option to be something else.".format(f))
        sys.exit()

# if the rename is in the args, build a dict
rename_dict = {}
if "rename" in args:
    for line in args.rename:
        line = line.strip()
        if line:
            fields = line.split("\t")
            rename_dict[fields[0]] = fields[1]

# get all the protein lengths
prot_to_len = {}
for record in SeqIO.parse(args.proteins, "fasta"):
    prot_to_len[record.id] = len(record.seq)

# get all the scaffold lengths with seqIO
scaf_to_len = {}
for record in SeqIO.parse(args.assembly, "fasta"):
    thisrecord = record.id
    if record.id in rename_dict:
        thisrecord = rename_dict[thisrecord]
    scaf_to_len[thisrecord] = len(record.seq)


# save all the protein IDs and their locations
prot_dicts = {}
# now get the protein ID lines
for line in args.gff:
    line = line.strip()
    if line and not line.startswith("#"):
        fields = line.split("\t")
        if (fields[2] == "CDS") and ("protein_id=" in fields[8]):
            protid = [x for x in fields[8].split(";")
                      if "protein_id=" in x][0].replace("protein_id=","")
            scaf = fields[0]
            if scaf in rename_dict:
                scaf = rename_dict[scaf]
            if scaf in scaf_to_len:
                direc = fields[6]
                start = int(fields[3])
                stop  = int(fields[4])
                if (protid not in prot_dicts) and (protid in prot_to_len):
                    prot_dicts[protid] = {
                        "scaf": scaf,
                        "direc": direc,
                        "start": start,
                        "stop": stop}
                else:
                    if start < prot_dicts[protid]["start"]:
                        prot_dicts[protid]["start"] = start
                    if stop > prot_dicts[protid]["stop"]:
                        prot_dicts[protid]["stop"] = stop

print("Number of proteins ", len(prot_to_len))
print("Length of prot_dicts ", len(prot_dicts))
# print a chrom file for the proteins
with open(chrom2loc, "w") as chromfile:
    df = pd.DataFrame.from_dict(prot_dicts ).T
    df.index.name = "protein"
    df = df.reset_index()
    df = df[["protein", "scaf", "direc", "start", "stop"]]
    # sort the df by scaf start stop
    df = df.sort_values(["scaf", "start", "stop"]).reset_index(drop=True)
    #now print the df with no headers to chromfile
    df.to_csv(chromfile, sep="\t", header=False, index=False)

# make a new file of proteins that appear in the genome
with open(proteins, "w") as protfile:
    args.proteins.seek(0)
    # parse through the proteins from the args and filter with prot_dicts
    for record in SeqIO.parse(args.proteins, "fasta"):
        if record.id in prot_dicts:
            # write with SeqIO
            SeqIO.write(record, protfile, "fasta")

# just copy the file args.assembly to assembly - don't do any special filtering
with open(assembly, "w") as outfile:
    args.assembly.seek(0)
    for line in args.assembly:
        print(line, file=outfile, end="")
