[![DOI](https://zenodo.org/badge/266213821.svg)](https://zenodo.org/badge/latestdoi/266213821)
## <a name="started"></a>Getting Started

```sh
# INSTALLATION
git clone https://github.com/conchoecia/genome_assembly_pipelines
cd genome_assembly_pipelines && make

# EXAMPLE RUN 
# We want to insert small scaffolds into chromosomes with Hi-C data.
#  First we must copy the example config file to the directory in which we want to run the program.
cp genome_assembly_pipelines/example_configs/config_GAP_sort_scaffolds_by_hic.yaml ./config.yaml
# then we run snakemake using the associated Snakefile. All the files are saved in the current directory.
snakemake --snakefile genome_assembly_pipelines/scripts/GAP_sort_scaffolds_by_hic_insert
```


## Table of Contents

- [Getting Started](#started)
- [Users' Guide](#uguide)
  - [Installation](#install)
  - [General usage](#general)
  - [Use cases](#cases)
    - [Insert small scaffolds into chromosome-scale scaffolds](#insert-scaffolds-hic)
  - [Getting help](#help)
  - [Citing GAP](#cite)

## <a name="uguide"></a>Users' Guide

This is a collection of `snakemake` pipelines that are designed to aid
in the _de novo_ assembly of genomes.

**Movitvation:** As much as we would like genome assembly to be a
one-command operation, it involves many steps, parameters to optimize,
and comparisons between assembly versions.  These tasks become more
complicated when working with highly heterozygous species.

Currently these scripts are largely undocumeted and they are 100%
unpublished. Documentation will be added below as needed for
collaborations or publications.

### <a name="install"></a>Installation

GAP was implemented in a linux environment with a bash shell, python
3, and snakemake. Every pipeline is implemented in snakemake, and
additional required programs are installed locally when executing
`make` or at runtime.

The only steps for installation are to clone the repository and make
the included programs.

```sh
git clone https://github.com/conchoecia/genome_assembly_pipelines
cd genome_assembly_pipelines && make
```

### <a name="general"></a>General usage

The general workflow for GAP requires the user to copy an example
configuration file into the directory in which they wish the program
to be run. I recommend doing this in an empty directory to keep the
analysis organized and to avoid overwriting and existing `config.yaml`
file.

```sh
cp genome_assembly_pipelines/example_configs/config_GAP_sort_scaffolds_by_hic.yaml ./config.yaml
```

You will find instructions for how to fill out the configuration
file inside the configuration file itself.

```sh
vim config.yaml
```

After setting up the configuration file, snakemake is used to run the
pipeline. You should refer to the snakemake documentation to determine
the number of cores to run, or how to configure snakemake for a
cluster environment like SLURM. Here, we add the options `-r -p` to
have more informative messages be printed to the terminal.

```sh
snakemake -r -p --snakefile genome_assembly_pipelines/scripts/GAP_sort_scaffolds_by_hic_insert
```

Output files are typically contained within a single folder of the
same name as the tool that was run.

### <a name="cases"></a>Use cases

These scripts all pertain to various operations used in genome
assembly.

#### <a name="insert-scaffolds-hic"></a>Insert small scaffolds into chromosome-scale scaffolds

If you have a chromosome-scale genome assembly with many small
scaffolds, it is desirable to insert those small scaffolds into the
chromosomes if there is some evidence to do so. Hi-C reads provide
such evidence.

This script maps Hi-C data to the genome assembly, and uses a rolling
window to determine the location in the chromsomes to which the
scaffold has the strongest Hi-C interaction signal. The user selects
what quantile of Hi-C interaction strength, and what minimum scaffold
size they wish to insert. The program applies these cutoffs as a logical
AND.

Insertions into the chromosomes are made only by adding the small
scaffodls between existing contigs. No contigs are broken in this
procedure. The input `.fasta` file and the output `.fasta` file have
the same number of contigs.

After the insertion process, there now remains small scaffolds that
are not inserted into the chromosome-scale scaffolds, but still have
Hi-C connections to the chromosome-scale scaffolds. These are sorted
in the fasta file based on their strongest Hi-C location in the
chromosome. These are all appended to the `.fasta` file after the
chromosome-scale scaffolds.

Finally, there may be scaffolds that have no Hi-C connections to the
chromosome-scale scaffolds. These are appended to the end of the
`.fasta` file in no particular order.

*Consideration 1*: If you run this program before running
[purge_dups][purge_dups] or [purge_haplotigs][purge_haplotigs], then
it is important to manually curate the genome assembly and remove
duplicate sequences based on the Hi-C map. If there is community
demand for a more explicit explanation, please add an issue in the
[github issues][issuepage] tab above.

*Consideration 2*: This program does not replace manual genome
curation, it only makes this step easier by removing the intial steps
of inserting scaffolds.

### <a name="help"></a>Getting help

This is currently the only documentation for these scripts. For
problems that arise please use the [github issues][issuepage] tab above.

### <a name="cite"></a>Citing GAP

There currently is no publication for citing the GAP package. Please
consider citing the github repository directly.

[issuepage]: https://github.com/conchoecia/genome_assembly_pipelines/issues
[purge_dups]: https://github.com/dfguan/purge_dups
[purge_haplotigs]: https://bitbucket.org/mroachawri/purge_haplotigs
