# `GAP_filter_competitive_map_coverage` — Quick Start Guide

## Inputs

You need three things in your `config.yaml`:

1. **Assembly** — a FASTA file that contains both chromosome-scale scaffolds and unplaced contigs
2. **Chromosome names** — the exact scaffold names from your FASTA that are chromosome-scale
3. **Long reads** — HiFi and/or ONT reads in FASTQ format

```yaml
assemblies:
  mygenome: "/path/to/assembly.fasta"

chromosomes:
  - Chr01
  - Chr02
  - Chr03

long_reads:
  hifi:
    file: "/path/to/hifi_reads.fastq.gz"
    type: "map-hifi"
```

The assembly does not need to be haplotype-resolved. You just need chromosome-scale scaffolds and unplaced contigs in the same FASTA file.

## Output

The main output is an **interactive HTML report** showing each unplaced contig plotted by:
- **X-axis:** % of the contig that aligns to chromosome scaffolds
- **Y-axis:** long-read coverage

You also get a CSV summary table with per-contig statistics.

### How to interpret the plot

- **High mapping % + high coverage:** likely haplotigs — redundant with chromosomes
- **Low mapping % + low coverage:** likely contaminants or artifacts
- **Low mapping % + high coverage:** potentially real unplaced sequences — investigate further
- **High mapping % + low coverage:** may be assembly artifacts or collapsed regions

Click on any point in the HTML report to see alignment dot plots showing exactly where the contig maps to chromosomes.

## Filtering the genome

First run **without filtering** to look at the HTML report and decide on thresholds:

```yaml
filter_redundant: False
trim_redundant_regions: False
```

Then enable filtering to produce a cleaned assembly:

```yaml
# Remove contigs where ≥95% aligns to chromosomes with MAPQ ≥30
filter_redundant: True
filter_min_coverage: 95.0
filter_min_mapq: 30

# Trim redundant regions from partially overlapping contigs
trim_redundant_regions: True
trim_min_block_size: 10000
trim_min_mapq: 20
trim_min_keep_piece: 5000
```

The final cleaned assembly is the `*_combined.fasta` file — chromosomes plus retained/trimmed contigs.

For full documentation, see [GAP_filter_competitive_map_coverage_README.md](GAP_filter_competitive_map_coverage_README.md).
