# Resource Estimates for GAP_hic_map7 Pipeline

This document summarizes the resource requirements (memory and runtime) for each rule in the GAP_hic_map7 snakemake pipeline.

## Resource Units
- **mem_mb**: Memory in megabytes (1000 MB = 1 GB)
- **runtime**: Time in minutes

## Rules and Their Resource Requirements

### Assembly Processing Rules

| Rule | Memory (MB) | Runtime (min) | Notes |
|------|------------|---------------|-------|
| `chrom_size` | 2,000 | 5 | Quick parsing of fasta headers |
| `compile_chromap` | 4,000 | 15 | One-time compilation |
| `check_assembly` | 4,000 | 10 | Validates unique sequence names |
| `index_ref` | Dynamic: max(8000, fasta_size/1000000 * 2) | 30 | Scales with assembly size |
| `generate_assembly_for_hic_gen` | 2,000 | 5 | Simple awk processing |
| `editable_assembly_file` | 4,000 | 10 | Python script processing |
| `genome_bed` | 2,000 | 5 | Quick bed file generation |
| `gaps_from_assembly` | 4,000 | 15 | Python parsing for N-gaps |

### Hi-C Mapping Rules

| Rule | Memory (MB) | Runtime (min) | Notes |
|------|------------|---------------|-------|
| `hic_to_pairs` | Dynamic: max(32000, fasta_size/1000000 * 8) | Dynamic: max(180, total_fastq_GB * 2) | Main Hi-C mapping step - most resource intensive |
| `gzip_pairs` | 4,000 | 30 | Compression of pairs file |
| `pairs2hiclongformat` | 8,000 | 30 | Format conversion |
| `index_pairs` | 4,000 | 15 | Pairix indexing |

### 3D-DNA and Juicebox Rules

| Rule | Memory (MB) | Runtime (min) | Notes |
|------|------------|---------------|-------|
| `get3ddna` | 2,000 | 10 | Git clone operation |
| `JBAT_pairs_to_hic` | 64,000 | 120 | Memory-intensive hic file generation |

### Matrix Generation and Processing Rules

| Rule | Memory (MB) | Runtime (min) | Notes |
|------|------------|---------------|-------|
| `generate_coverage_stats_final_pairs` | 8,000 | 30 | Python script analyzing pairs |
| `make_bins_individual` | 16,000 | 60 | Cooler matrix generation |
| `diagnostic_plot` | 8,000 | 20 | Plot generation |
| `normalize_matrix` | 12,000 | 30 | Matrix normalization |
| `balance_matrix` | 16,000 | 45 | Cooler balancing operation |
| `zoomify_normalized_and_balanced_matrix` | 12,000 | 30 | Multi-resolution aggregation |
| `zoomify_matrix` | 12,000 | 30 | Multi-resolution aggregation |
| `make_hic_matrix_simple` | 12,000 | 30 | Format conversion |
| `best_connected_scaffold` | 16,000 | 45 | Python analysis of connections |

### Pretext Rules

| Rule | Memory (MB) | Runtime (min) | Notes |
|------|------------|---------------|-------|
| `make_pretextmap` | Dynamic: max(32000, fasta_size/1000000 * 8) | Dynamic: max(180, total_fastq_GB * 2) | Similar to hic_to_pairs |

### Gap Analysis Rules

| Rule | Memory (MB) | Runtime (min) | Notes |
|------|------------|---------------|-------|
| `make_gaps_beddb` | 8,000 | 20 | HiGlass database creation |

### Telomere and K-mer Analysis Rules

| Rule | Memory (MB) | Runtime (min) | Notes |
|------|------------|---------------|-------|
| `telomere_kmer_positions_bigwig` | 4,000 | 20 | K-mer position tracking |
| `telomere_kmer_bigwig_to_multivec` | 6,000 | 15 | Multi-track file creation |
| `kmer_positions_bigwig` | 4,000 | 20 | K-mer position tracking |
| `telomere_lookup_longreads` | 4,000 | 30 | Grep-based search in reads |
| `filter_telo_bams_longreads` | 8,000 | 30 | Picard filtering |
| `filter_telodir_mapdir` | 4,000 | 15 | Samtools filtering |
| `bedgraph_of_telomeres` | 8,000 | 30 | Coverage calculation |
| `bigwig_of_telomeres` | 4,000 | 10 | Format conversion |
| `telomere_bigwig_to_multivec` | 6,000 | 15 | Multi-track file creation |
| `bedgraph_and_bed2bigwig` | 4,000 | 10 | Generic format conversion |

### Transcript Mapping Rules

| Rule | Memory (MB) | Runtime (min) | Notes |
|------|------------|---------------|-------|
| `symlink_the_transcripts` | 1,000 | 2 | Symlink creation |
| `map_transcripts_to_genome` | Dynamic: max(16000, fasta_size/1000000 * 4) | Dynamic: max(60, transcript_GB * 30) | Minimap2 splice alignment |
| `index_transcripts_bams` | 2,000 | 10 | Samtools indexing |
| `get_coverage_of_transcripts_bams` | 8,000 | 30 | Coverage calculation |

### Long Read Mapping Rules

| Rule | Memory (MB) | Runtime (min) | Notes |
|------|------------|---------------|-------|
| `symlink_the_LRs` | 1,000 | 2 | Symlink creation |
| `map_LR_to_genome` | Dynamic: max(16000, fasta_size/1000000 * 4) | Dynamic: max(60, LR_GB * 30) | Minimap2 long read alignment |
| `index_LR_bams` | 2,000 | 10 | Samtools indexing |
| `get_coverage_of_LR_bams` | 8,000 | 30 | Coverage calculation |

## Dynamic Resource Calculations

Several rules use dynamic resource allocation based on input file sizes:

### Memory Scaling
1. **Assembly-based scaling**: `max(base_memory, assembly_size_MB * multiplier)`
   - Used for indexing and mapping operations
   - Ensures sufficient memory for large genomes

2. **Read-based scaling**: Calculated from fastq.gz file sizes
   - Hi-C mapping scales with total read data
   - Long read mapping scales with individual file sizes

### Runtime Scaling
1. **Read mapping**: `max(base_time, file_size_GB * time_per_GB)`
   - Hi-C: ~2 minutes per GB of compressed reads
   - Long reads: ~30 minutes per GB of compressed reads
   - Transcripts: ~30 minutes per GB of compressed reads

2. **Fixed time operations**: Most format conversions and simple analyses

## Recommendations

### For Small Genomes (< 500 MB)
- Most rules will use base memory allocations
- Hi-C mapping: ~32-40 GB RAM, 3-6 hours
- Total pipeline: 4-12 hours with adequate cores

### For Medium Genomes (500 MB - 2 GB)
- Hi-C mapping: ~40-80 GB RAM, 6-12 hours
- JBAT conversion: 64 GB RAM required
- Total pipeline: 12-24 hours with adequate cores

### For Large Genomes (> 2 GB)
- Hi-C mapping: 80+ GB RAM, 12+ hours
- Consider splitting Hi-C libraries for parallel processing
- JBAT may require 64-128 GB RAM
- Total pipeline: 24-48 hours with adequate cores

## Notes on Resource Estimation

1. **Dynamic resources** are calculated at runtime based on input file sizes
2. **Thread counts** are set per rule and may affect memory usage
3. **Temporary files** are automatically cleaned but require disk space during execution
4. **The bottleneck** is typically the `hic_to_pairs` rule for Hi-C mapping
5. **Peak memory usage** occurs during `JBAT_pairs_to_hic` (64 GB minimum)

## Slurm Integration

These resources can be used directly with Snakemake's cluster execution:

```bash
snakemake --cluster "sbatch --mem={resources.mem_mb} --time={resources.runtime}" \
          --default-resources mem_mb=4000 runtime=30
```

Or in a cluster profile configuration file.
