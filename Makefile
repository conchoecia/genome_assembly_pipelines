CPPF=-g -O2 -W -Wall -Wno-unused-parameter -Wpointer-arith -Wshadow -Wundef -std=c++11
FSD=src/fstats/

all: bin bin/fasta_stats dependencies/SALSA/run_pipeline.py dependencies/SALSA/break_contigs_start bin/pilon-1.23.jar

bin:
	mkdir bin

bin/fasta_stats: ${FSD}fasta_stats.cc ${FSD}itoa.cc ${FSD}open_compressed.h
	g++ ${CPPF} -o bin/fasta_stats ${FSD}*.cc

bin/pilon-1.23.jar: bin
	wget --directory-prefix bin 'https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar'

dependencies/SALSA/run_pipeline.py:
	mkdir -p dependencies
	cd dependencies; git clone https://github.com/marbl/SALSA.git
	for file in dependencies/SALSA/*.py; do sed -i '1s;^;#!/usr/bin/python2\n;' $${file}; done

dependencies/SALSA/break_contigs_start: dependencies/SALSA/run_pipeline.py
	cd dependencies/SALSA; make
