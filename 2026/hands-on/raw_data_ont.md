# Nanopore Basecalling and Methylation Analysis Pipeline

## Overview

Perform Oxford Nanopore basecalling with three accuracy models (SUP, HAC, FAST), quality assessment, and modified base calling with methylation analysis.

## Prerequisites

### Software Requirements

- **Dorado** v1.3.1: Oxford Nanopore basecaller (`/project/varcall_training/bin/dorado-1.3.1-linux-x64/bin/dorado`)
- **Modkit**: Modified base analysis tool (`/project/varcall_training/bin/varcall_latest.sif`)
- **R** scripts for plotting (`/project/varcall_training/data/real/src/qual.r`, `/project/varcall_training/data/real/src/meth.r`)
- **AWK**: For binning quality scores (`/project/varcall_training/data/real/src/bin.awk`)

### Input Files

- **POD5 files**: Raw nanopore data (`/project/varcall_training/data/real/sample/sample.pod5`)
- **Reference genome**: FASTA format (`/project/varcall_training/data/real/reference/chr11.fa`)
- **Basecalling models**:
  - `/project/varcall_training/data/real/models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0` (super-accurate)
  - `/project/varcall_training/data/real/models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0` (high-accurate)
  - `/project/varcall_training/data/real/models/dna_r10.4.1_e8.2_400bps_fast@v4.3.0` (fast)
- **Modified base model**: `/project/varcall_training/data/real/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v2.0.1/`

## Pipeline Steps

### 0. Make your own directory

```bash
cd ~
mkdir -p varcall_training/day1/ont
cd varcall_training/day1/ont
#start an interactive session with GPU
```

### 1. Standard Basecalling (SUP Model)

```bash
/project/varcall_training/bin/dorado-1.3.1-linux-x64/bin/dorado basecaller \
  -x cuda:all \
  -o results/sup_results \
  /project/varcall_training/data/real/models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
  /project/varcall_training/data/real/sample/sample.pod5
```

**Parameters**:
- `-x cuda:all`: Use all available GPUs
- `-o`: Output directory for BAM files
- Model: Super-accurate basecalling model
- Input: POD5 files

**Output**: Unmapped BAM (`uBAM`)

### 2. Generate Sequencing Summary (SUP)

```bash
/project/varcall_training/bin/dorado-1.3.1-linux-x64/bin/dorado summary \
  results/sup_results/FFPE_v2_v_PrePCR_AP_29042023/PrePCR/20230429_1600_2E_PAO83395_124388f5/bam_pass/PAO83395_pass_124388f5_0a73e955_0.bam \
  > results/sup_results/sequencing_summary.tsv
```

**Output**: Tab-separated file containing read-level statistics (quality scores, length, etc.)

### 3. Bin Quality Scores (SUP)

```bash
tail -n+2 results/sup_results/*tsv | \
  cut -f 16 | \
  awk -v delta='0.5' -f /project/varcall_training/data/real/src/bin.awk \
  > results/sup_results/qual.tsv
```

**Parameters**:
- `tail -n+2`: Skip header
- `cut -f 16`: Extract quality score column
- `delta='0.5'`: Bin size for quality score histogram

**Output format (example)**:
```text
0.0 0.5 0
0.5 1.0 0
...
17.5 18.0 6
18.0 18.5 10
```

### 4. HAC Model Basecalling

```bash
/project/varcall_training/bin/dorado-1.3.1-linux-x64/bin/dorado basecaller \
  -x cuda:all \
  -o results/hac_results \
  /project/varcall_training/data/real/models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0 \
  /project/varcall_training/data/real/sample/sample.pod5

/project/varcall_training/bin/dorado-1.3.1-linux-x64/bin/dorado summary \
  results/hac_results/FFPE_v2_v_PrePCR_AP_29042023/PrePCR/20230429_1600_2E_PAO83395_124388f5/bam_pass/PAO83395_pass_124388f5_0a73e955_0.bam \
  > results/hac_results/sequencing_summary.tsv

tail -n+2 results/hac_results/*tsv | \
  cut -f 16 | \
  awk -v delta='0.5' -f /project/varcall_training/data/real/src/bin.awk \
  > results/hac_results/qual.tsv
```

### 5. FAST Model Basecalling

```bash
/project/varcall_training/bin/dorado-1.3.1-linux-x64/bin/dorado basecaller \
  -x cuda:all \
  -o results/fast_results \
  /project/varcall_training/data/real/models/dna_r10.4.1_e8.2_400bps_fast@v4.3.0 \
  /project/varcall_training/data/real/sample/sample.pod5

/project/varcall_training/bin/dorado-1.3.1-linux-x64/bin/dorado summary \
  results/fast_results/FFPE_v2_v_PrePCR_AP_29042023/PrePCR/20230429_1600_2E_PAO83395_124388f5/bam_pass/PAO83395_pass_124388f5_0a73e955_0.bam \
  > results/fast_results/sequencing_summary.tsv

tail -n+2 results/fast_results/*tsv | \
  cut -f 16 | \
  awk -v delta='0.5' -f /project/varcall_training/data/real/src/bin.awk \
  > results/fast_results/qual.tsv
```

### 6. Quality Score Comparison Visualization

```bash
Rscript /project/varcall_training/data/real/src/qual.r \
  results/sup_results/qual.tsv \
  results/hac_results/qual.tsv \
  results/fast_results/qual.tsv \
  results
```

**Outputs**:
- `results/quality_comparison_faceted.png`: Side-by-side faceted comparison

### 7. Modified Base Calling with Reference Alignment

```bash
/project/varcall_training/bin/dorado-1.3.1-linux-x64/bin/dorado basecaller \
  -x cuda:all \
  -o results/sup_mod_results \
  /project/varcall_training/data/real/models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
  --modified-bases-models /project/varcall_training/data/real/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v2.0.1/ \
  /project/varcall_training/data/real/sample/sample.pod5 \
  --reference /project/varcall_training/data/real/reference/chr11.fa
```

**Parameters**:
- `--modified-bases-models`: Model for detecting 5mC and 5hmC modifications
- `--reference`: Reference genome for alignment (anchors reads to genomic positions)

**Output**: BAM file with `MM` and `ML` tags containing modification information

### 8. Generate bedMethyl File with Modkit

```bash
singularity exec -B "$PWD,/project/varcall_training,/localscratch" /project/varcall_training/bin/varcall_latest.sif modkit pileup \
  results/sup_mod_results/FFPE_v2_v_PrePCR_AP_29042023/PrePCR/20230429_1600_2E_PAO83395_124388f5/bam_pass/PAO83395_pass_124388f5_0a73e955_0.bam \
  results/sup_mod_results/meth.bed \
  --ref /project/varcall_training/data/real/reference/chr11.fa \
  --combine-strands \
  --cpg \
  --modified-bases 5mC 5hmC
```

**Parameters**:
- `--ref`: Reference genome
- `--cpg`: Report only CpG dinucleotide positions
- `--combine-strands`: Combina calls on forward and reverse strand
- `--log-filepath`: Output log file for QC

**Output**: bedMethyl file with 18 columns containing modification counts and fractions

### 9. Visualize Methylation Profile

Run:

```bash
Rscript /project/varcall_training/data/real/src/meth.r results/sup_mod_results/meth.bed
```

**Output** - `results/sup_mod_results/meth_methylation_profile.png`: methylation profile of the region for the sequenced sample

## Output File Formats

### bedMethyl Format
The modkit output contains 18 columns:

| Column | Name | Description |
|--------|------|-------------|
| 1 | chrom | Chromosome name |
| 2 | start | 0-based start position |
| 3 | end | 0-based exclusive end |
| 4 | mod_code | Modification code (m=5mC, h=5hmC) |
| 5 | score | Valid coverage (N_valid_cov) |
| 6 | strand | + or - |
| 7-9 | start/end/color | BED compatibility fields |
| 10 | N_valid_cov | Total valid calls |
| 11 | frac_modified | Modification fraction (%) |
| 12 | N_mod | Modified base calls |
| 13 | N_canonical | Canonical base calls |
| 14 | N_other_mod | Other modification calls |
| 15 | N_delete | Deletion count |
| 16 | N_fail | Failed calls |
| 17 | N_diff | Different base calls |
| 18 | N_nocall | No-call count |


## Directory Structure

```
results/
|-- fast_results
|   |-- FFPE_v2_v_PrePCR_AP_29042023
|   |   `-- PrePCR
|   |       `-- 20230429_1600_2E_PAO83395_124388f5
|   |           `-- bam_pass
|   |               `-- PAO83395_pass_124388f5_0a73e955_0.bam
|   |-- qual.tsv
|   `-- sequencing_summary.tsv
|-- hac_results
|   |-- FFPE_v2_v_PrePCR_AP_29042023
|   |   `-- PrePCR
|   |       `-- 20230429_1600_2E_PAO83395_124388f5
|   |           `-- bam_pass
|   |               `-- PAO83395_pass_124388f5_0a73e955_0.bam
|   |-- qual.tsv
|   `-- sequencing_summary.tsv
|-- quality_comparison_faceted.png
|-- sup_mod_results
|   |-- FFPE_v2_v_PrePCR_AP_29042023
|   |   `-- PrePCR
|   |       `-- 20230429_1600_2E_PAO83395_124388f5
|   |           `-- bam_pass
|   |               |-- PAO83395_pass_124388f5_0a73e955_0.bam
|   |               `-- PAO83395_pass_124388f5_0a73e955_0.bam.bai
|   |-- meth.bed
|   `-- meth_methylation_profile.png
`-- sup_results
    |-- FFPE_v2_v_PrePCR_AP_29042023
    |   `-- PrePCR
    |       `-- 20230429_1600_2E_PAO83395_124388f5
    |           `-- bam_pass
    |               `-- PAO83395_pass_124388f5_0a73e955_0.bam
    |-- qual.tsv
    `-- sequencing_summary.tsv

20 directories, 14 files
```

## References

- [Dorado Documentation](https://github.com/nanoporetech/dorado)
- [Modkit Documentation](https://github.com/nanoporetech/modkit)
- [bedMethyl Format Specification](https://www.encodeproject.org/data-standards/wgbs/)
