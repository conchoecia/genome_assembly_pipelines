#!/usr/bin/env python3

"""
author: darrin t. schultz
github: @conchoecia
date:   2023-12-21

This file converts an assembly file from artisanal to fasta format.
The assembly file is one that shows how the contigs are turned into scaffolds.

Removes terminal Ns from the sequences.
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys

def parse_args():
    """
    The args are:
      -a --assembly_file - The path to the assembly file. We must check if it exists.
      -f --fasta_file    - The path to the fasta file. We must check if it exists.
      -g --gap_length    - The length of the gap between the contigs. Default is 100.
      -n --names         - The names for the first N scaffolds of the assembly. Use a space-separated list.
      -p --prefix        - The prefix for the scaffold names. Default is "scaffold".
      -h --help          - Print out all of the options
    """
    parser = argparse.ArgumentParser(description="Converts an assembly file to a fasta file.")
    parser.add_argument("-a", "--assembly_file", required = True, type=str, help="The path to the assembly file.")
    parser.add_argument("-f", "--fasta_file",    required = True, type=str, help="The path to the fasta file.")
    parser.add_argument("-g", "--gap_length",    required = False, type=int, help="The length of the gap between the contigs. Default is 100.", default=100)
    parser.add_argument("-n", "--first_n_names",         required = False, type=str, help="The names for the first N scaffolds of the assembly. Use a space-separated list.")
    parser.add_argument("-p", "--prefix",        required = False, type=str, help="The prefix for the scaffold names. Default is 'scaffold'.", default="scaffold")
    args = parser.parse_args()

    # check if the assembly file and the fasta file exist
    if os.path.isfile(args.assembly_file) == False:
        print("The assembly file does not exist.")
        exit()
    if os.path.isfile(args.fasta_file) == False:
        print("The fasta file does not exist.")
        exit()
    # if the names are given, parse the string into a list of strings
    if args.first_n_names:
        args.first_n_names = args.first_n_names.split(" ")
    else:
        args.first_n_names = []

    return args

# set up a class called AssemblyFileFasta that will be manipulated to create the fasta file.
class AssemblyFileFasta:
    """
    This class will be used to create the fasta file from the assembly file.

    The assembly file is split into two parts. The first part is a description of how the contigs are turned into scaffolds.
    The second part of the assembly file is the way in which the fragments above must be oriented to create the scaffolds.

    >scaffold_1:::fragment_1 1 137844079
    >scaffold_1:::fragment_2 2 153358929
    >scaffold_2:::fragment_1 3 228461863
    >scaffold_3:::fragment_1 4 176027019
    .
    .
    .
    -2 -1
    -7
    -12 -11 -10
    16 17
    """
    def __init__(self, assembly_file, fasta_file):
        self.assembly_file = assembly_file
        self.fasta_file    = fasta_file
        # This dictionary contains the instructions on how to break up the fasta file into fragments
        self.scaffold_to_fragments = self.parse_assembly_file_to_coordinates()
        # This list of lists contains the instructions on how to generate the final genome assembly that will be printed
        self.assembly_instructions = self.parse_assembly_file_to_assembly_instructions()
        # This contains the entire fasta file as fragments. The whole thing will be loaded into memory. Use a machine with >10GB of RAM.
        self.fasta_fragments = {} # update directly from the function to avoid making a copy of the genome in RAM
        self.parse_fasta_file_into_fragments()

    def parse_assembly_file_to_coordinates(self):
        """
        This reads the assembly file and parses it into an easy-to-use dictionary that will be used later.

        This only parses this part of the file:
        >scaffold_1:::fragment_1 1 137844079
        >scaffold_1:::fragment_2 2 153358929
        >scaffold_2:::fragment_1 3 228461863
        """
        assemdict = {}
        handle = open(self.assembly_file, "r")
        for line in handle:
            line = line.strip()
            if line: # some lines may be blank
                if line[0] == ">" and ":::" in line:
                    fields    = line.split(" ")
                    scaf      = fields[0].split(":::")[0].replace(">","")
                    fragment  = fields[0].split(":::")[1]
                    fragindex = int(fields[1])
                    seqlen    = int(fields[2])
                    # check if the scaffold is in the dictionary
                    if scaf not in assemdict:
                        assemdict[scaf] = {}
                    assemdict[scaf][fragindex] = seqlen
        handle.close()
        return assemdict

    def parse_assembly_file_to_assembly_instructions(self):
        """
        This reads in the assembly file and parses the instructions on how to take the fragments and create scaffolds.
        The sequences will be printed in the same order as they appear in the assembly file.

        This only parses this part of the assembly file:
        -2 -1
        -7
        -12 -11 -10

        The output is a list of lists that looks like this:
        [[-2, -1], [-7], [-12, -11, -10]]
        """
        assembly_instructions = []
        handle = open(self.assembly_file, "r")
        for line in handle:
            line = line.strip()
            if line: # some lines may be blank
                # only get the other lines
                if line[0] not in [">","#"] and ":::" not in line:
                    fields = [int(x) for x in line.split(" ")]
                    assembly_instructions.append(fields)
        handle.close()
        return assembly_instructions

    def parse_fasta_file_into_fragments(self):
        """
        This reads in a fasta file and converts it into a dictionary of fragments.
        The indices are the fragment numbers from column 2 of the assembly file.
        The values are BioPython Seq objects.
        """
        # parse through each sequence in the fasta file
        for seq in SeqIO.parse(self.fasta_file, "fasta"):
            # check if the scaffold is in the dictionary
            if seq.id not in self.scaffold_to_fragments:
                raise IOError("We have seen a scaffold in the fasta file that was not in the assembly file. Assembly in fasta file: {}".format(seq.id))
            # get the fragments for this scaffold, go through them in ascending numeric order
            prev_position = 0 # we use prev_position and curr_position to keep track of where we are in the sequence
            curr_position = 0
            for fragindex in sorted(self.scaffold_to_fragments[seq.id].keys()):
                curr_position += self.scaffold_to_fragments[seq.id][fragindex]
                fragseq = seq[prev_position:curr_position]
                # check if the fragment is in the dictionary
                if fragindex not in self.fasta_fragments:
                    self.fasta_fragments[fragindex] = fragseq
                else:
                    raise IOError("We have seen a fragment in the fasta file that was already in the dictionary. Fragment: {}".format(fragindex))
                prev_position = curr_position

    def _trim_terminal_Ns(self, seq):
        """
        This trims terminal Ns or ns from a sequence.
        """
        return seq.strip("N").strip("n")

    def generate_fasta(self, gapsize, outprefix, first_n_names):
        """
        Takes all of the information that has been read so far and generates a fasta file.
        """
        # make a string for the gap
        gapstring = "".join(["N"] * gapsize)

        # Here, we keep track of the first N scaffold names.
        # This will be whittled down as we print out scaffolds.
        namebank = [x for x in first_n_names]

        # go through each scaffold in the assembly instructions
        # The counter is only used for the scaffold counter after the first N names have been consumed.
        counter = 1
        for thisscaf in self.assembly_instructions:
            # Temoporarily put the raw strings into this list.
            # They will be joined together with gaps of Ns later
            new_seqs = []
            # go through each fragment in the scaffold
            for thisfrag in thisscaf:
                # check if the fragment is in the dictionary
                if abs(thisfrag) not in self.fasta_fragments:
                    raise IOError("We have seen a fragment in the assembly instructions that was not in the fasta file. Fragment: {}".format(thisfrag))
                # If the fragment is negative, we need to reverse complement it, otherwise, just add it to the list.
                # We only add the string value.
                if thisfrag < 0:
                    new_seqs.append(self._trim_terminal_Ns(str(self.fasta_fragments[abs(thisfrag)].reverse_complement().seq)))
                else:
                    new_seqs.append(self._trim_terminal_Ns(str(self.fasta_fragments[thisfrag].seq)))
            # join the sequences together with Ns
            outsequence = gapstring.join(new_seqs)
            # create a new SeqRecord object with the new outsequence, and for the newid the outprefix and the counter
            newid = ""
            if len(namebank) > 0:
                newid = namebank.pop(0)
            else:
                newid = outprefix + "_" + str(counter)
                counter += 1
            seqrec = SeqIO.SeqRecord(Seq(outsequence), id=newid, description="")
            # print the sequence to sys.stderr
            SeqIO.write(seqrec, sys.stdout, "fasta")

def main():
    args = parse_args()
    assembly = AssemblyFileFasta(args.assembly_file, args.fasta_file)
    # generate the fasta file
    assembly.generate_fasta(args.gap_length, args.prefix, args.first_n_names)

if __name__ == "__main__":
    main()