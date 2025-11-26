# GAP_phase Multi-Library Support - Implementation Summary

## Overview
The GAP_phase pipeline has been updated to support multiple sequencing libraries per data type (Hi-C, long reads, Chicago, Illumina). 

**KEY CHANGE: All libraries now use pre-aligned BAM files instead of FASTQ files.**

This allows you to:
1. Combine multiple Hi-C libraries from different restriction enzymes (DpnII, MboI, Arima, etc.)
2. Merge multiple PacBio or ONT sequencing runs
3. Process multiple Chicago or Illumina libraries together
4. **Reuse BAM files** across multiple phasing runs without re-mapping
5. **Each library is processed through sambamba markdup -r to remove duplicates**

## Key Changes

### 1. Configuration Format - ALL BAM FILES
**New format** (recommended):
```yaml
libraries:
  hic:
    DpnII:
      bam: "/path/to/DpnII_aligned.bam"
    MboI:
      bam: "/path/to/MboI_aligned.bam"
  
  longreads:
    PacBio_run1:
      bam: "/path/to/pacbio1.bam"
    PacBio_run2:
      bam: "/path/to/pacbio2.bam"
```

**Old format** (still supported for backward compatibility):
```yaml
HIC_BAM: "/path/to/hic.bam"
LONGREAD_BAM: "/path/to/longread.bam"
CHI_BAM: "/path/to/chicago.bam"
ILL_BAM: "/path/to/illumina.bam"
```

### 2. Pipeline Workflow

The new workflow is:

```
For each library:
  1. Create softlinks/symlinks to input BAM files

For each data type (hic, chi, ill, longreads):
  2. Merge all libraries of that type → Single merged BAM
  3. Index the merged BAM
  4. Split by chromosome
  5. Add read groups
  6. Remove duplicates with sambamba markdup -r
  7. Continue with HapCUT2 phasing...
```

### 3. How to Prepare Your BAM Files

**IMPORTANT for Hi-C and Chicago reads:**
Hi-C BAM files MUST be aligned with BWA using `-5`, `-S`, `-P` flags as recommended by HapCUT2:
```bash
bwa mem -t 16 -5SPM reference.fasta hic_R1.fq.gz hic_R2.fq.gz | \
  samtools view -hb | samtools sort > hic_aligned.bam
```

Reference: https://github.com/vibansal/HapCUT2/issues/122

**For Chicago reads** (similar to Hi-C):
```bash
bwa mem -t 16 -5SPM reference.fasta chi_R1.fq.gz chi_R2.fq.gz | \
  samtools view -hb | samtools sort > chicago_aligned.bam
```

**For Illumina reads** (standard alignment):
```bash
bwa mem -t 16 reference.fasta ill_R1.fq.gz ill_R2.fq.gz | \
  samtools view -hb | samtools sort > illumina_aligned.bam
```

**For PacBio reads:**
```bash
minimap2 -t 16 -ax map-pb reference.fasta pacbio.fq.gz | \
  samtools view -hb | samtools sort > pacbio_aligned.bam
```

**For Oxford Nanopore reads:**
```bash
minimap2 -t 16 -ax map-ont reference.fasta ont.fq.gz | \
  samtools view -hb | samtools sort > ont_aligned.bam
```

### 4. New Rules

- `make_library_bam_softlink`: Creates symlinks for BAM files per library
- `merge_libraries`: Merges multiple libraries of the same type into one BAM
- **Removed**: All FASTQ alignment rules (no longer needed)

### 5. Duplicate Removal

The `mark_duplicates_in_split_bams` rule uses `sambamba markdup -r`:

```bash
# Removes duplicates (not just marks them):
sambamba markdup -r -t 4 input.bam output.bam
```

Benefits:
- **Removes duplicates** instead of just marking them (cleaner data for phasing)
- Faster execution (multithreaded)
- Lower memory footprint
- No metrics file needed

### 6. File Structure

```
GAP_phase/
├── input/
│   ├── assembly/
│   │   └── input.fasta
│   └── bams/
│       ├── hic.DpnII.bam
│       ├── hic.MboI.bam
│       ├── longreads.PacBio_run1.bam
│       └── longreads.PacBio_run2.bam
├── output/
│   └── bams/
│       └── temp/                  # Merged BAMs (one per type)
│           ├── hic_sorted.bam
│           ├── longreads_sorted.bam
│           └── ...
```

## Example Configurations

See the example config files:
- `example_configs/config_GAP_phase_multi_library.yaml` - Full example with multiple libraries
- `example_configs/config_GAP_phase_OLD_vs_NEW.yaml` - Comparison of old vs new formats

## Backward Compatibility

The old configuration format is still supported. The pipeline will automatically detect which format you're using and convert old configs to the new internal format.

## Usage Examples

### Example 1: Multiple Hi-C Libraries
```yaml
libraries:
  hic:
    DpnII_rep1:
      bam: "/data/hic_dpnii_1_aligned.bam"
    DpnII_rep2:
      bam: "/data/hic_dpnii_2_aligned.bam"
    Arima:
      bam: "/data/hic_arima_aligned.bam"
```

### Example 2: Multiple PacBio Runs
```yaml
libraries:
  longreads:
    SMRT_cell_1:
      bam: "/data/pacbio_cell1.bam"
    SMRT_cell_2:
      bam: "/data/pacbio_cell2.bam"
    SMRT_cell_3:
      bam: "/data/pacbio_cell3.bam"
```

### Example 3: Mixed Data Types
```yaml
libraries:
  hic:
    Arima:
      bam: "/data/arima_aligned.bam"
  
  longreads:
    PacBio:
      bam: "/data/pacbio.bam"
    ONT:
      bam: "/data/nanopore.bam"
  
  ill:
    NovaSeq:
      bam: "/data/illumina_aligned.bam"
```

## Benefits of BAM-Only Approach

1. **Reusability**: Align reads once, reuse BAMs for multiple phasing runs with different parameters
2. **Speed**: Skip time-consuming alignment step when re-running pipeline
3. **Consistency**: Use the same alignment parameters across all your analyses
4. **Flexibility**: Mix and match different sequencing runs easily
5. **Storage**: Keep compressed BAM files instead of larger FASTQ files

## Testing

To test the new format with your existing data:

1. Align your reads to BAM files using the commands shown above
2. Create your config file with BAM paths (see examples)
3. Run: `snakemake -s snakefiles/GAP_phase --configfile config.yaml -n`
4. Check that all expected outputs are listed
5. Run the pipeline

## Notes

- Library names must be unique within each type
- Library names can contain letters, numbers, and underscores
- If a library type is not needed, simply omit it from the config
- All BAMs must be coordinate-sorted and aligned to the same reference assembly
- sambamba must be installed (it's now in the Makefile: `make bin/sambamba`)
- **Hi-C and Chicago BAMs must be aligned with `-5SPM` flags for HapCUT2 compatibility**
