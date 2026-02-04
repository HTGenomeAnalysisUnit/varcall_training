# Alignment and Assembly: Short vs. Long Reads

This practical exercise demonstrates the **fundamental differences between short-read and long-read sequencing** for genomic analysis.

The focus spans two complementary analyses:

1. **Reference-based alignment**: How reads map to known sequence
2. **De novo assembly**: Reconstructing sequence without a reference

Participants will:
- Align GIAB short- and long-read data to reference
- Assemble both read types independently 
- Compare assemblies to the reference to assess contiguity and accuracy
- Visualize differences in complex genomic regions

## 1. Conceptual Background

### 1.1 Short Reads vs. Long Reads

| Property | Short Reads (Illumina) | Long Reads (PacBio HiFi) |
|----------|------------------------|--------------------------|
| Length | 100-300 bp | 10-25 kb |
| Error rate | <0.1% | <1% (HiFi: <0.5%) |
| Throughput | Very high | Moderate |
| Cost per base | Low | Higher |
| Repetitive regions | Difficult | Can span repeats |
| Assembly contiguity | Fragmented (low N50) | Highly contiguous |

### 1.2 Why Assembly Matters

**Reference-based alignment** works well for:

- Variant calling in well-characterized regions
- Known gene analysis
- Population studies

**De novo assembly** is essential for:

- Structural variants (insertions, deletions, duplications)
- Novel sequences not in reference
- Resolving paralogous gene families
- Understanding haplotype structure

### 1.3 Assembly Quality Metrics

- **N50**: The contig length at which 50% of the assembly is in contigs of that size or larger (higher = better contiguity)
- **Number of contigs**: Fewer contigs = more contiguous assembly
- **Assembly length**: Total bases assembled (should approximate reference length)
- **Mapping rate**: Percentage of assembly that aligns to reference
- **Identity**: Sequence accuracy when aligned to reference

## 2. Training Regions

Create a BED file with regions demonstrating different assembly challenges:

```bash
cat > training_regions.bed << 'EOF'
chr17   72062001   76562000   control_region
chr22    42031864  42245566  CYP2D_locus
chr4   143685957   144248249   GYP_locus
EOF
```

**Region characteristics**:
- `control_unique_region`: Unique, low-complexity - both read types assemble well
- `*_locus`: Locus with known copy-number variable events


## 3. Directory Setup

Navigate to a folder for which you have read/write permission:

```bash
mkdir -p reference short-reads long-reads
mkdir -p alignment_ex/{aligned,assembled,analysis}
```

## 4. Data Preparation

### 4.1 Subset GIAB Data to Training Regions

**Long reads** (PacBio HiFi):

```bash
samtools view -@ 8 -b -M -L training_regions.bed \
  -o long-reads/HG002.subset.bam \
  --write-index \
  https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam
```

**Short reads** (Illumina 2x250 bp):

```bash
samtools view -@ 8 -b -M -L training_regions.bed \
  -o short-reads/HG002.subset.bam \
  --write-index \
  https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam
```

**Flags**:
- `-M`: Keep mate pairs even if one falls outside target, use index to subset alignment to regions
- `-L`: Restrict to BED regions
- `--write-index`: Create CSI index simultaneously

### 4.2 Extract Chr20 Reference

```bash
cat /project/varcall_training/data/partial/genome/hg38/chr{4,17,22}.fa > reference/reference.fa
samtools faidx reference/reference.fa
```

### 4.3 Convert BAM to FASTQ

**Short reads** (paired-end):

```bash
samtools collate -@ 8 -O short-reads/HG002.subset.bam short-reads/tmp_collated | \
  samtools fastq -@ 8 -N \
    -1 short-reads/hg002.r1.fq.gz \
    -2 short-reads/hg002.r2.fq.gz \
    -s short-reads/hg002.sng.fq.gz \
    -0 /dev/null
```

**Long reads** (single-end):

```bash
samtools fastq -@ 8 -0 long-reads/hg002.fq.gz \
  long-reads/HG002.subset.bam
```

### 4.4 Read Statistics

```bash
# Short reads
echo "Short read pairs:"
zcat short-reads/hg002.r1.fq.gz | echo $((`wc -l`/4))

echo "Short read N50:"
seqkit stats -T -N 50 short-reads/hg002.r1.fq.gz

# Long reads  
echo "Long reads:"
zcat long-reads/hg002.fq.gz | echo $((`wc -l`/4))

echo "Long read N50:"
seqkit stats -T -N 50 long-reads/hg002.fq.gz 
```

## 5. Part 1: Reference-Based Alignment

### 5.1 Build Minimap2 Index

```bash
minimap2 -d reference/reference.mmi reference/reference.fa
```

### 5.2 Align Short Reads (Default Parameters)

```bash
cd alignment_ex/aligned

minimap2 \
  -a -x sr -t 8 \
  -R "@RG\tID:HG002_SR\tPL:ILLUMINA\tSM:HG002" \
  ../../reference/reference.mmi \
  ../../short-reads/hg002.r1.fq.gz \
  ../../short-reads/hg002.r2.fq.gz | \
  samtools sort -@ 8 --write-index -o hg002.sr.bam
```

**Preset `-x sr`**: k=21, w=11, optimized for Illumina 100-300bp reads

### 5.3 Align Long Reads (Default Parameters)

```bash
minimap2 \
  -a -x map-hifi -t 8 \
  -R "@RG\tID:HG002_HIFI\tPL:PACBIO\tSM:HG002" \
  ../../reference/reference.mmi \
  ../../long-reads/hg002.fq.gz | \
  samtools sort -@ 8 --write-index -o hg002.hifi.bam
```

**Preset `-x map-hifi`**: k=19, w=10, optimized for PacBio HiFi reads

### 5.4 Alignment Statistics

```bash
cd ../analysis

for bam in ../aligned/*.bam; do
  echo "=== $(basename "$bam") ==="
  samtools flagstat "$bam"
  echo
done > alignment_stats.txt

cat alignment_stats.txt
```

### 5.5 Mapping Quality in each region Region

```bash
BED=../../training_regions.bed
SR_BAM=../aligned/hg002.sr.bam
LR_BAM=../aligned/hg002.hifi.bam

echo "=== Short reads: MAPQ summary per region ==="
while read -r chr start end label; do
  [[ "$chr" =~ ^# ]] && continue
  region="${chr}:${start}-${end}"
  samtools view "$SR_BAM" "$region" | \
    awk -v lbl="$label" -v reg="$region" '
      {mapq[$5]++; total++}
      END {
        if (total == 0) {
          printf "%-30s %-25s  Total reads: %6d  MAPQ60: 0.0%%  MAPQ<60: 0.0%%\n", lbl, reg, total;
        } else {
          m60 = mapq[60]+0;
          below60 = total - m60;
          pct60 = (m60/total)*100;
          pctb = (below60/total)*100;
          printf "%-30s %-25s  Total reads: %6d  MAPQ60: %5.1f%%  MAPQ<60: %5.1f%%\n", lbl, reg, total, pct60, pctb;
        }
      }'
done < "$BED"

echo
echo "=== Long reads: MAPQ summary per region ==="
while read -r chr start end label; do
  [[ "$chr" =~ ^# ]] && continue
  region="${chr}:${start}-${end}"
  samtools view "$LR_BAM" "$region" | \
    awk -v lbl="$label" -v reg="$region" '
      {mapq[$5]++; total++}
      END {
        if (total == 0) {
          printf "%-30s %-25s  Total reads: %6d  MAPQ60: 0.0%%  MAPQ<60: 0.0%%\n", lbl, reg, total;
        } else {
          m60 = mapq[60]+0;
          below60 = total - m60;
          pct60 = (m60/total)*100;
          pctb = (below60/total)*100;
          printf "%-30s %-25s  Total reads: %6d  MAPQ60: %5.1f%%  MAPQ<60: %5.1f%%\n", lbl, reg, total, pct60, pctb;
        }
      }'
done < "$BED"
```

**Expected observation**: Short reads will show more MAPQ=0, while long reads have higher MAPQ due to spanning unique flanking sequence.

## 6. Part 2: De Novo Assembly

### 6.1 Assemble Short Reads with SPAdes

SPAdes is optimized for short-read assembly using de Bruijn graphs.

```bash
cd ../assembled

spades.py \
  --isolate \
  -1 ../../short-reads/hg002.r1.fq.gz \
  -2 ../../short-reads/hg002.r2.fq.gz \
  -s ../../short-reads/hg002.sng.fq.gz \
  -o spades_assembly \
  -t 8 -m 16

# Copy final contigs
cp spades_assembly/contigs.fasta hg002.spades.contigs.fa
```

**Key flags**:
- `--isolate`: Optimized for high-coverage isolate sequencing (not metagenomic)
- `-t 8`: Use 8 threads
- `-m 16`: Use up to 16 GB RAM

**Note**: SPAdes will try multiple k-mer sizes internally (automatic k-mer selection).

### 6.2 Assemble Long Reads with Hifiasm

Hifiasm is designed for PacBio HiFi reads and produces phased, haplotype-resolved assemblies.

```bash
hifiasm \
  -o hg002.hifiasm \
  -t 8 \
  ../../long-reads/hg002.fq.gz

# Convert GFA to FASTA (primary contigs)
awk '/^S/{print ">"$2"\n"$3}' hg002.hifiasm.bp.p_ctg.gfa > hg002.hifiasm.contigs.fa
```

**Output files**:
- `.bp.p_ctg.gfa`: Primary contigs (one haplotype representation)
- `.bp.hap1.p_ctg.gfa` / `.bp.hap2.p_ctg.gfa`: Haplotype-resolved assemblies (if enough heterozygosity)

### 6.3 Assembly Statistics

```bash
cd ../analysis

# Basic stats
for asm in ../assembled/*.contigs.fa; do
  echo "=== $(basename "$asm") ==="
  seqkit stats -T -N 50 "$asm" | column -t
done > assembly_stats.txt

cat assembly_stats.txt
```

**Compare**:
- **Number of contigs**: SPAdes will have many more (fragmented assembly)
- **N50**: Hifiasm will have much higher N50 (more contiguous)
- **Max length**: Hifiasm contigs can span 100+ kb

## 7. Assembly-to-Reference Comparison

### 7.1 Align Assemblies to Reference

**Short-read assembly**:

```bash
cd ../aligned

minimap2 -ax asm5 -t 8 \
  ../../reference/reference.mmi \
  ../assembled/hg002.spades.contigs.fa | \
  samtools sort -@ 8 --write-index -o hg002.spades.vs_ref.bam
```

**Long-read assembly**:

```bash
minimap2 -ax asm5 -t 8 \
  ../../reference/reference.mmi \
  ../assembled/hg002.hifiasm.contigs.fa | \
  samtools sort -@ 8 --write-index -o hg002.hifiasm.vs_ref.bam
```

**Preset `-x asm5`**: For assembly-to-reference alignment (~5% divergence)

### 7.2 Assembly Coverage and Identity

```bash
cd ../analysis

for bam in ../aligned/*.vs_ref.bam; do
  echo "=== $(basename "$bam") ==="
  samtools coverage "$bam"
  echo
done > assembly_coverage.txt

cat assembly_coverage.txt
```

**Metrics**:
- `coverage`: Fraction of reference covered by assembly
- `meandepth`: How many assembly contigs overlap each reference position
- `meanmapq`: Average mapping quality of assembly contigs

### 7.3 Assembly Contiguity by Region

Extract contigs mapping to each training region:

```bash
#!/bin/bash

BED=../../training_regions.bed

> assembly_by_region.txt

while read -r chr start end label; do
  [[ "$chr" =~ ^# ]] && continue
  region="${chr}:${start}-${end}"
  
  echo "=== Region: $label ($region) ===" >> assembly_by_region.txt
  
  echo "Short-read assembly (SPAdes):" >> assembly_by_region.txt
  samtools view -c ../aligned/hg002.spades.vs_ref.bam "$region" >> assembly_by_region.txt
  
  echo "Long-read assembly (Hifiasm):" >> assembly_by_region.txt
  samtools view -c ../aligned/hg002.hifiasm.vs_ref.bam "$region" >> assembly_by_region.txt
  
  echo "" >> assembly_by_region.txt
done < "$BED"

cat assembly_by_region.txt
```

**Interpretation**: SPAdes will have many small contigs; Hifiasm may have longer contigs

### 8 Dotplot Comparison

Visualize assembly-to-reference alignment structure:

```bash
# Requires minimap2 + R/Python plotting
minimap2 -cx asm5 \
  ../../reference/reference.mmi \
  ../assembled/hg002.hifiasm.contigs.fa > hifiasm.paf

minimap2 -cx asm5 \
  ../../reference/reference.mmi \
  ../assembled/hg002.spades.contigs.fa > spades.paf
```
There are packages one can use to visualize the generated `paf` files as dotplots (`pafR` in the R environment, for instance)

## 9. Software Requirements

All the tools are available in the dedicated singularity container (`/project/varcall_training/bin/varcall_latest.sif`)
