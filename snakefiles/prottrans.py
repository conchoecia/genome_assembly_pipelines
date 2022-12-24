#! /usr/bin/env python
#
# v3.3 rearranged code, removed junk 2015-03-18
# v3.2 fixed general output number bug 2015-02-20
# v3.1 fixed 6-frame output number bug 2015-01-27
# v3.0 removed revcomp, some major reworks 2014-04-18
# v2.2 functionalized statistic read out for orf measurements 2014-04-17
# v2.1 added no-rename mode, commented a bit, fixed multiframe bug 2014-02-06
# v2.0 remade the functions, cleaned up code, indexing functions 2013-11-09
# v1.8 added specification of frame by user 04/feb/2013
# v1.7 added alternate translation tables
# v1.6 added methionine flag to enable m-searching in all modes 12/jul/2012
# v1.5 forked with revcomp, default translate mode is 'm' 25/may/2012
# v1.4 corrected minor formatting errors 2/apr/2012
# v1.3 converted to argparse, added functions 9/feb/2012
# v1.2 added 6-frame translation 7/feb/2012
# v1.1 added translation tool 3/feb/2012
# v1.0 simple reverse complementer of a fasta file 4/jan/2012

'''prottrans.py last modified 2015-10-05
    a protein translator requiring Bio python

  for fast protein translation, choosing the best protein (most cases)
prottrans.py sequences.fasta > proteins.fasta

  to print all six frames
prottrans.py -t 6 sequences.fasta > 6frame_translated.fasta
  to instead keep the full frame with the longest ORF
prottrans.py -t f sequences.fasta > fulltranslated.fasta
  to print the nucleotides with the translation
prottrans.py -t t sequences.fasta > nucl_and_translation.txt

  ## for options related to orf length counting

  to print out lengths of all orfs and not translate
prottrans.py -t o sequences.fasta > orfcountss.txt
  to print out stats of all max-length-normalized orfs
prottrans.py -t l sequences.fasta > l-norm-orfs.txt
  to print out stats of all transcript-length-normalized orfs
prottrans.py -t p sequences.fasta > t-norm-orfs.txt

  to show the list of genetic codes and quit
prottrans.py -c 0 -
'''

import sys
import time
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from collections import deque

codelist='''List of Genetic Codes

1    The Standard Code
2    The Vertebrate Mitochondrial Code
3    The Yeast Mitochondrial Code
4    The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
5    The Invertebrate Mitochondrial Code
6    The Ciliate, Dasycladacean and Hexamita Nuclear Code
9    The Echinoderm and Flatworm Mitochondrial Code
10   The Euplotid Nuclear Code
11   The Bacterial, Archaeal and Plant Plastid Code
12   The Alternative Yeast Nuclear Code
13   The Ascidian Mitochondrial Code
14   The Alternative Flatworm Mitochondrial Code
15   Blepharisma Nuclear Code
16   Chlorophycean Mitochondrial Code
21   Trematode Mitochondrial Code
22   Scenedesmus Obliquus Mitochondrial Code
23   Thraustochytrium Mitochondrial Code
24   Pterobranchia Mitochondrial Code'''

def get_three_frame(seq_obj,code):
	allprotframes = []
	intronphase = deque([None,-2,-1])
	# some reason there is a problem subsetting 0:0 if entire seq should be translated
	intronphase.rotate(len(seq_obj) % 3)
	allprotframes.append(seq_obj[0:intronphase[0]].translate(table=code))
	allprotframes.append(seq_obj[1:intronphase[1]].translate(table=code))
	allprotframes.append(seq_obj[2:intronphase[2]].translate(table=code))
	return allprotframes

def get_six_frame(seq_record,code, noreverse):
	'''does the 6-frame translation and returns a list of all proteins from six frames'''
	dna_seq = seq_record.seq
	allprotframes = get_three_frame(dna_seq,code)
	if not noreverse:
		allprotframes.extend( get_three_frame(dna_seq.reverse_complement(),code) )
	return allprotframes

def get_met_orfs(preorfs):
	'''list, if M is present, get orf from position of the first M to the end of orf'''
	orfs = [orf[orf.find("M"):] for orf in preorfs if orf.count('M')]
	return orfs

def sizecut_orfs(preorfs, above, below):
	orfs = [orf for orf in preorfs if below > len(orf) >= above]
	return orfs

def frame_to_orfs(protframe):
	orfs = protframe.split('*')
	return orfs

def get_longest_orf(protframe, requiremet, above, below):
	'''returns a tuple of length of the longest protein, longest protein sequence
    and integer count of the fragments that are greater than 0'''
	orfs = frame_to_orfs(protframe)
	orfnum = len(orfs)
	if requiremet:
		orfs = get_met_orfs(orfs)
	orfs = sizecut_orfs(orfs, above, below)
	# sort the orfs by length
	orfs.sort(key=len, reverse=True)
	fragmentlengths = [len(f) for f in orfs]
	# count the number of orfs of length zero
	realfrags = orfnum - fragmentlengths.count(0)
	# return the tuple
	try:
		return (fragmentlengths[0], orfs[0], realfrags)
	# this is to account for an empty list, as calling fragmentlengths[0]
	except IndexError:
	# as a complicated bit of code with max() was used before, this was changed from ValueError
		return (0, '', 0)

def get_orf_stats(everyorf):
	everyorf.sort()
	sumoforfs = sum(everyorf)
	# remove orfs of length zero from average and median
	orfavg = float(sumoforfs)/(len(everyorf)-everyorf.count(0))
	orfmed = everyorf[(len(everyorf)//2)+everyorf.count(0)]
	return orfavg, orfmed

def full_frame_outprot(transmode, sixframeprots, keptprot):
	'''keptprot is a tuple of the longest prot and the frame
    so for transmode f keep the full frame of that index
    otherwise keep the protein itself'''
	if transmode == 'f':
		outprot = sixframeprots[keptprot[1]] # this is the index of the best prot
	else: # implying if transmode == 'n':
		outprot = keptprot[0] 
	return outprot

def prep_outprotein(outseqrec, protein, outid):
	'''this function only updates values for the SeqRecord'''
	outseqrec.seq = protein
	outseqrec.id = outid
	# due to parsing errors, the name is printed twice if id and description do not match
	# so changed to blank string
	outseqrec.description = ""
	# void function

def get_frame_id(sequenceid, i):
	outidstring = "{}_{}".format(sequenceid,i)
	return outidstring

def check_frame(frame):
	if abs(frame) > 3 or frame==0:
		print >> sys.stderr, "WARNING: FRAME {} IS NOT VALID, MUST USE +1,2,3, or -1,-2,-3".format(frame)
		print >> sys.stderr, "EXITING", time.asctime()
		sys.exit()
	else:
		return frame

def get_python_frame(inframe):
	if inframe > 0:
		# get the python index from frames 1 to 3, giving 0 to 2
		outframe = inframe-1
	else:
		# otherwise as frames are from -1 to -3, abs()+2 gives 3 to 5
		outframe = abs(inframe)+2
	return outframe

def get_best_prot(sixframeprots, requiremet, above, below):
	'''a list for each of length of the longest protein, longest protein sequence, 
    and integer count of the fragments that are greater than 0'''
	fragmentmax = []
	longestprots = []
	# for each frame in the sequences, get the longest orf
	for protframe in sixframeprots:
		longestbyframe = get_longest_orf(protframe, requiremet, above, below)
		fragmentmax.append(longestbyframe[0])
		longestprots.append(longestbyframe[1])
	overallmax = max(fragmentmax)
	if overallmax:
		bestprots = [(x[1],i) for i,x in enumerate(zip(fragmentmax,longestprots)) if x[0]==overallmax]
	else:
		bestprots = []
	# returns a list of tuples containing the sequence and the frame number
	return bestprots

def derive_normalization(transmode,protframe):
	if transmode == 'l':
		# normalize by the length of the longest ORF
		normalizefactor = float(max([len(o) for o in frame_to_orfs(protframe)]))
	elif transmode == 'p':
		# normalize by the total length of the frame
		normalizefactor = float(len(str(protframe)))
	else: # else meaning all letters other than l or p
		# print the raw length, hence divided by 1
		normalizefactor = 1
	return normalizefactor

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', type = argparse.FileType('rU'), default = '-', help="fasta format file, or - for stdin")
	parser.add_argument('-a', '--above', type=float, metavar='N', default=1.0, help="only count sequences longer than N")
	parser.add_argument('-b', '--below', type=float, metavar='N', default=1000000000.0, help="only count sequences shorter than N")
	parser.add_argument('-c', '--code', type=int, default=1, help="use alternate genetic code, type 0 to print code list and exit")
	parser.add_argument('-f', '--frame', type=int, help="get translation frame of 1,2,3,-1,-2,-3")
	parser.add_argument('-t', '--translate', choices='nf6tolp', help="type of translation, default: n", default='n')
	parser.add_argument('-m', '--methionine', action='store_false', help="do not require a methionine")
	parser.add_argument('-n', '--norename', action='store_true', help="forbid renaming of sequences, discard multiframe")
	parser.add_argument('-r', '--no-reverse', action='store_true', help="only translate forward strand")
	parser.add_argument('-s', '--statistics', action='store_true', help="do additional calculations")
	parser.add_argument('-v', '--verbose', action='store_true', help="display more output")
	args = parser.parse_args(argv)

	if args.code==0:
		sys.exit(codelist)

	transmode = args.translate
	norename = args.norename

	print >> sys.stderr, "# Translating sequences from {}".format(args.input_file.name), time.asctime()
	if transmode == 'f' or transmode == '6' or transmode == 'n':
		do_sequences = True
		if args.no_reverse:
			print >> sys.stderr, "# Translating only forward strand", time.asctime()
	else:
		print >> sys.stderr, "# Calculating stats for {}".format(args.input_file.name), time.asctime()
		do_sequences = False
		everyorf = [0] # why is this a list of zero?

	# basic counts
	aacountsum = 0
	readcount = 0
	writecount = 0
	nuclsum = 0
	nullcount = 0
	multiframes = 0

	if args.frame:
		useframe = get_python_frame(check_frame(args.frame))

	startclock=time.asctime()
	starttime= time.time()

	for seq_record in SeqIO.parse(args.input_file, "fasta"):
		readcount += 1
		nuclsum += len(seq_record.seq)
		# in all cases, start with 6-frame translation
		sixframeprots = get_six_frame(seq_record,args.code, args.no_reverse)
		if do_sequences:
			# must create base id string, otherwise numbers are appended each frame
			baseid = str(seq_record.id)
			# if printing all 6 frames, then ignore all future processing
			if transmode == '6':
				multiframes += 1
				for i, frame in enumerate(sixframeprots):
					writecount += 1
					aacountsum += len(frame)
					outprotid = get_frame_id(baseid,i)
					prep_outprotein(seq_record, frame, outprotid)
					wayout.write(seq_record.format("fasta"))
			# otherwise start checking for other features
			else:
				# useframe only has value if it passed the check
				if args.frame:
					keepframes = [(sixframeprots[useframe], useframe)]
				else:
					keepframes = get_best_prot(sixframeprots, args.methionine, args.above, args.below)
				# test if there are any good frames
				if len(keepframes):
					if len(keepframes) > 1:
						if args.verbose:
							print >> sys.stderr, "%d ORFS IN: %s" % (len(keepframes),outprotid)
						# multiframes counts the number of prots with multiple frames
						# instead of the total proteins derived from this, which would be 2 or more
						multiframes += 1
					for keptprot in keepframes:
						# for cases of duplicates or two short ORFs
						if norename:
							outprotid = baseid
							# if renaming is turned off, then disregard multiframe proteins
							if len(keepframes) > 1:
								if args.verbose:
									print >> sys.stderr, "IGNORING REDUNDANT PROTEIN {}".format(outprotid)
								# end for loop here and continue with the next protein
								continue
						else:							
							outprotid = get_frame_id(baseid,keptprot[1])
						# determines to use the entire frame or not
						outprot = full_frame_outprot(transmode, sixframeprots, keptprot)
						protlen = len(outprot)
						# updates information for the seq record
						prep_outprotein(seq_record, outprot, outprotid)
						wayout.write(seq_record.format("fasta"))
						writecount += 1
						aacountsum += protlen
						if args.verbose:
							if len(keepframes) > 1:
								print >> sys.stderr, "MULTIFRAME PROTEIN: {}AAs: {}".format(protlen, outprotid)
							else:
								print >> sys.stderr, "FOUND PROTEIN {}".format(outprotid)
				# for cases where no prots were kept, probably because there was no M or too short
				else:
					# all valid prots have been length filtered, or no M and stop in any frame
					if args.verbose:
						print >> sys.stderr, "NO VALID PROTS FOUND: {}".format(seq_record.id)
					# invalid - no sequences were found that contained both M and stop in any frame
					nullcount += 1
		else:
			# for printing length of all orfs of all frames
			for protframe in sixframeprots:
				# using alternate normalization modes
				normalizefactor = derive_normalization(transmode,protframe)
				allorfs = frame_to_orfs(protframe)
				for orf in allorfs:
					writecount += 1
					orflen = len(str(orf))
					normorflength = orflen/normalizefactor
					print >> wayout, normorflength
					if args.statistics:
						everyorf.append(normorflength)

	print >> sys.stderr, "Counted %d records in %.2f minutes" % (readcount, (time.time()-starttime)/60)
	print >> sys.stderr, "Total counted bases: {}".format(nuclsum)
	print >> sys.stderr, "Printed %d sequences; %d are multiframe" % (writecount, multiframes)
	print >> sys.stderr, "Total translated aminoacids is: {}".format(aacountsum)
	print >> sys.stderr, "Number of invalid sequences: {}".format(nullcount)
	if args.statistics and not do_sequences:
		orfmean, orfmedian = get_orf_stats(everyorf)
		print >> sys.stderr, "Median is: {} and mean is: {:.2f}".format(orfmedian,orfmean)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)

