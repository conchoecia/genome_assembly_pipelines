#!/usr/bin/env python

"""
In this script, take the fasta file of protein sequences, and filter a gff file
  to only include the proteins that are in the fasta file. This script also generates
  a new chrom file that only includes this file subset.

The parsing works by looking in the ID= field of the gff file and matching it to
  the protein fasta headers.

The output chrom format is headerless: {protein ID, scaffold, strand, start, end}

XP_033732128.1  NC_047015.1     +       65488   122595
XP_033757380.1  NC_047015.1     +       171889  173437
XP_033732147.1  NC_047015.1     -       176365  227045
XP_033732158.1  NC_047015.1     +       227717  296214
XP_033732177.1  NC_047015.1     -       304704  342850
XP_033732169.1  NC_047015.1     -       304704  342850
XP_033757397.1  NC_047015.1     -       363415  363415
XP_033732191.1  NC_047015.1     -       404433  634172
"""

import argparse
import sys

def parse_args():
    """
    The args we need:
      - protein fasta file
      - gff file
      - outprefix (used to name the output chrom and gff files)
    """
    parser = argparse.ArgumentParser(description='Filter gff file to only include proteins in fasta file')
    parser.add_argument("-f", "--fasta", type=str, help='Protein fasta file')
    parser.add_argument("-g", "--gff", type=str, help='Gff file')
    parser.add_argument("-o", "--outprefix", type=str, help='Output prefix for chrom and gff files')
    args = parser.parse_args()
    # if no args, print help. This didn't work, try something else
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def read_fasta(fasta_file):
    """Read protein IDs from a FASTA file"""
    protein_ids = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                protein_id = line.split()[0][1:]  # Extract ID from header
                protein_ids.add(protein_id)
    return protein_ids

def filter_gff(gff_file, protein_ids, output_gff):
    """Filter GFF file to only include entries with matching protein IDs"""
    with open(gff_file, 'r') as infile, open(output_gff, 'w') as outfile:
        for line in infile:
            if line.startswith("##gff-version 3"):
                outfile.write(line)
                continue
            if line.startswith('#'):
                # we don't want these in the final gff file
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                print("Erroneous line: ", line)
                raise ValueError(f"Line in GFF file has fewer than 9 fields: {line}")
            attributes = fields[8]
            for attr in attributes.split(';'):
                if attr.startswith('ID=') or attr.startswith('Parent='):
                    protein_id = attr.split('=')[1]
                    if protein_id in protein_ids:
                        outfile.write(line)
                        break

def generate_chrom_file(filtered_gff, output_chrom):
    """Generate a chrom file from the filtered GFF file.
    These will be sorted based on scaffold, then by start, then by stop position
    """
    entries = []
    with open(filtered_gff, 'r') as infile, open(output_chrom, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                raise ValueError(f"Line in GFF file has less than 9 fields: {line}")
            scaffold = fields[0]
            feature_type = fields[2]
            start = fields[3]
            end = fields[4]
            strand = fields[6]
            attributes = fields[8]
            protein_id = None
            for attr in attributes.split(';'):
                if attr.startswith('ID='):
                    protein_id = attr.split('=')[1]
                    break
            if protein_id and feature_type == "mRNA":
                entries.append((protein_id, scaffold, strand, start, end))
        # Now that we have all the entries sort them
        entries.sort(key=lambda x: (x[1], int(x[3]), int(x[4])))
        # Write to the output file
        for entry in entries:
            outfile.write('\t'.join(entry) + '\n')

def main():
    args = parse_args()
    protein_ids = read_fasta(args.fasta)
    output_gff   = f"{args.outprefix}.gff"
    output_chrom = f"{args.outprefix}.chrom"

    filter_gff(args.gff, protein_ids, output_gff)
    print(f"Filtered GFF file saved to {output_gff}")

    generate_chrom_file(output_gff, output_chrom)
    print(f"Chrom file saved to {output_chrom}")

if __name__ == '__main__':
    main()
