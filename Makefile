CPPF=-g -O2 -W -Wall -Wno-unused-parameter -Wpointer-arith -Wshadow -Wundef -std=c++11
FSD=src/fstats/

all: bin/fasta_stats dependencies/SALSA/run_pipeline.py \
     dependencies/SALSA/break_contigs_start \
     bin/bbmap/reformat.sh \
     bin/pilon-1.23.jar dependencies/pairix/bin/pairix \
     bin/fqjt bin/minlen_pair \
	 dependencies/miniprot

bin/fqjt: src/fq-jt.c
	mkdir -p bin
	gcc -gdwarf-2 -g src/fq-jt.c -lz -o bin/fqjt

bin/bbmap/reformat.sh: dependencies/bbtools/BBMap_38.90.tar.gz
	tar -C bin/ -zxvf dependencies/bbtools/BBMap_38.90.tar.gz

bin/minlen_pair: src/minlen_pair.c src/kseq.h
	mkdir -p bin
	gcc src/minlen_pair.c -lz -o bin/minlen_pair

bin/fasta_stats: ${FSD}fasta_stats.cc ${FSD}itoa.cc ${FSD}open_compressed.h
	mkdir -p bin
	g++ ${CPPF} -o bin/fasta_stats ${FSD}*.cc

bin/pilon-1.23.jar:
	mkdir -p bin
	wget --directory-prefix bin 'https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar'

dependencies/SALSA/run_pipeline.py:
	mkdir -p dependencies
	cd dependencies; git clone https://github.com/marbl/SALSA.git
	for file in dependencies/SALSA/*.py; do sed -i '1s;^;#!/usr/bin/python2\n;' $${file}; done

dependencies/SALSA/break_contigs_start: dependencies/SALSA/run_pipeline.py
	cd dependencies/SALSA; make

dependencies/pairix/setup.py:
	mkdir -p dependencies
	cd dependencies; git clone https://github.com/4dn-dcic/pairix

dependencies/pairix/bin/pairix: dependencies/pairix/setup.py
	cd dependencies/pairix; make

dependencies/miniprot:
	mkdir -p dependencies
	cd dependencies; git clone https://github.com/lh3/miniprot.git
	cd dependencies/miniprot; make