CPPF=-g -O2 -W -Wall -Wno-unused-parameter -Wpointer-arith -Wshadow -Wundef -std=c++11
FSD=src/fstats/

all: bin/fasta_stats dependencies/SALSA/run_pipeline.py \
     dependencies/SALSA/break_contigs_start \
     bin/bbmap/reformat.sh \
     bin/pilon-1.23.jar dependencies/pairix/bin/pairix \
     bin/fqjt bin/minlen_pair \
     dependencies/miniprot bin/sambamba \
     bin/k8 dependencies/hickit/hickit.js bin/fastp

bin/fqjt: src/fq-jt.c
	mkdir -p bin
	gcc -gdwarf-2 -g src/fq-jt.c -lz -o bin/fqjt

bin/bbmap/reformat.sh: dependencies/bbtools/BBMap_38.90.tar.gz
	tar -C bin/ -zxvf dependencies/bbtools/BBMap_38.90.tar.gz
	touch bin/bbmap/reformat.sh

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

bin/sambamba:
	@echo "Installing sambamba for Linux..."
	@mkdir -p bin
	@SAMBAMBA_URL=$$(curl -s https://api.github.com/repos/biod/sambamba/releases/latest | grep "browser_download_url.*linux-amd64-static" | cut -d '"' -f 4); \
	echo "Downloading from: $$SAMBAMBA_URL"; \
	wget -O bin/sambamba.gz "$$SAMBAMBA_URL"; \
	gunzip -f bin/sambamba.gz; \
	chmod +x bin/sambamba
	@echo "sambamba installed successfully to bin/sambamba"

bin/k8:
	@echo "Installing k8 JavaScript shell..."
	@mkdir -p bin
	@echo "Downloading k8..."
	@rm -f k8-1.2.tar.bz2
	wget https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2
	@echo "Extracting k8..."
	tar -jxf k8-1.2.tar.bz2
	@echo "Copying k8 binary..."
	cp -f k8-1.2/k8-x86_64-Linux bin/k8
	chmod +x bin/k8
	@echo "Verifying installation..."
	@ls -lh bin/k8
	@file bin/k8
	@echo "k8 installed successfully to bin/k8"
	@echo "Note: k8-1.2/ directory left in place for hickit compilation"

dependencies/hickit/hickit.js: bin/k8
	@echo "Installing hickit (requires k8)..."
	@mkdir -p dependencies
	cd dependencies; git clone https://github.com/lh3/hickit.git
	cd dependencies/hickit; make
	@echo "hickit installed successfully to dependencies/hickit/"

bin/fastp:
	@echo "Installing fastp for Linux..."
	@mkdir -p bin
	@echo "Downloading latest fastp binary..."
	wget -O bin/fastp http://opengene.org/fastp/fastp
	chmod +x bin/fastp
	@echo "fastp installed successfully to bin/fastp"


