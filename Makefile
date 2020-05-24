CPPF=-g -O2 -W -Wall -Wno-unused-parameter -Wpointer-arith -Wshadow -Wundef -std=c++11
FSD=src/fstats/

all: bin bin/fasta_stats dependencies/salsa/run_pipeline.py

bin:
	mkdir bin

bin/fasta_stats: ${FSD}fasta_stats.cc ${FSD}itoa.cc ${FSD}open_compressed.h
	g++ ${CPPF} -o bin/fasta_stats ${FSD}*.cc

dependencies/salsa/run_pipeline.py:
	mkdir -p dependencies
	cd dependencies
	git clone https://github.com/marbl/SALSA.git
	cd ../
