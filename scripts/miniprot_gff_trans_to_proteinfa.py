#!/usr/bin/env python

"""
Program: miniprot_gff_trans_to_proteinfa.py

This program takes a gff file from miniprot, output with the trans option, and exports
 the protein sequences in fasta format. The header of the protein sequence is the MP number,
 followed by a space, then the original protein sequence.

For this output:

##gff-version 3
##PAF	XP_033724493.1	368	119	368	+	OZ024763.1	30684678	14746944	14747688	498	783	0	AS:i:765	ms:i:765	np:i:193	fs:i:0	st:i:0	da:i:939	do:i:0	cg:Z:124M1I1M1I6M1I15M1D13M1I4M1D5M1I4M2I4M6I18M6D11M4D31M	cs:Z:*aacD:3*gatY*ctgS:2*aacD:4*tatF:13*agcG*gagD:12*gtcI*aacS:6*aagR:33*tacF:2*gttL:3*tacF*gaaN*ggtD:13*aaaR:4*acgA*acgA:3*cacP:6*tccN+H:1+G*accP*ttcY:4+Q:3*aacP:4*tccP*tacP*gcgY*cacS:3-cac:2*ccaS:2*acgA:1*aacT*ccgH*tccN*tccA*aacS:1+H*gccQ:3-cac:3*gtgI:1+N*tacF*agcN:2+TQ*cagA:3+PPGSTL:1*accS*ttcC:1*gccY*tacN:1*acaL*gccP*gccS*gccG*atgV:1*tcgL:1*atgP:1*ggcA-atgcgtatcgcccagccg*ccgG:1*aacT:1*ccaS:5*tcgN-ccggcgaccgcc*gccS*acgV:2*atgS*tcgG:2*gccP:1*accS:1*cagY*caaR:1*gccQ*gagG*ccgE:1*ttcI:2*atgV*cccH:4*gcaG:2
##STA	NSNKDLDPNQKPPYSYVALIAMAIKESSEKRLTLSGIYQYIVNKFPYYEKNKKGWQNSIRHNLSLNECFVKVPREGGGERKGNYWTVDPAYEGMFEKGNYRRRRRMKRPYRTTISLHKPLFADSCTFNQFALTKNYFSPSYAHQYSHQYPPWTLNPSSNAAGMGHMSQVSYSSCQRVPSTFGAYPTAAAMQSSMSGMRIAQPPTNYPQLNDYSPATAATSPMSPFAFTPQQQAEPSFNAMPYTYWAER*
OZ024763.1	miniprot	mRNA	14746945	14747691	765	+	.	ID=MP000001;Rank=1;Identity=0.6360;Positive=0.7395;Target=XP_033724493.1 120 368
OZ024763.1	miniprot	CDS	14746945	14747691	765	+	0	Parent=MP000001;Rank=1;Identity=0.6360;Target=XP_033724493.1 120 368
OZ024763.1	miniprot	stop_codon	14747689	14747691	0	+	0Parent=MP000001;Rank=1

The resulting fasta entry will be:
>MP000001 XP_033724493.1
NSNKDLDPNQKPPYSYVALIAMAIKESSEKRLTLSGIYQYIVNKFPYYEKNKKGWQNSIR
HNLSLNECFVKVPREGGGERKGNYWTVDPAYEGMFEKGNYRRRRRMKRPYRTTISLHKPL
FADSCTFNQFALTKNYFSPSYAHQYSHQYPPWTLNPSSNAAGMGHMSQVSYSSCQRVPST
FGAYPTAAAMQSSMSGMRIAQPPTNYPQLNDYSPATAATSPMSPFAFTPQQQAEPSFNAM
PYTYWAER*
"""

import sys

def parse_gff(file_path):
    sequences = {}
    current_protein = None
    current_id = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('##PAF'):
                parts = line.split('\t')
                if len(parts) > 1:
                    current_protein = parts[1]
                continue

            if line.startswith('##STA'):
                sequence = line.split('\t')[1] if '\t' in line else line[6:]
                if current_id and current_protein:
                    sequences[current_id] = (current_protein, sequence)
                continue

            fields = line.split('\t')
            if len(fields) < 9:
                continue

            attributes = {attr.split('=')[0]: attr.split('=')[1] for attr in fields[8].split(';') if '=' in attr}

            if 'ID' in attributes:
                current_id = attributes['ID']

    return sequences

def write_fasta(sequences, output_file):
    with open(output_file, 'w') as out:
        for mp_id, (protein, sequence) in sequences.items():
            out.write(f">{mp_id} {protein}\n")
            for i in range(0, len(sequence), 60):
                out.write(sequence[i:i+60] + "\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: miniprot_gff_trans_to_proteinfa.py <input.gff> <output.fasta>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    sequences = parse_gff(input_file)
    write_fasta(sequences, output_file)

if __name__ == '__main__':
    main()
