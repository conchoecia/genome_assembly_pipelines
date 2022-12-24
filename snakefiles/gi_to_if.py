#!/usr/bin/env python

import sys

chrom_to_count = {}

counter = 1
for line in sys.stdin:
    line = line.strip()
    if line:
        splitd = line.split("\t")
        col1 = "\t".join(splitd[0:3])
        col2 = "\t".join(splitd[3:6])
        count = splitd[-1]
        if col2 not in chrom_to_count:
            chrom_to_count[col2] = counter
            counter += 1
        else:
            pass
        print("{}\t{}\t{}".format(chrom_to_count[col1], chrom_to_count[col2], count))
