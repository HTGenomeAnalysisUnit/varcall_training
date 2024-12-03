# Alignment

In this chapter we start from raw FASTQ files from short and long read sequencing experiment and we perform alignment of the reads to the reference using [minimap2](https://github.com/lh3/minimap2). Sam-to-bam conversion is performed through [samtools](https://github.com/samtools/samtools).

## Input data

### A

Data folder: `/project/varcall_training/data/partial/fastq/illumina/`
Files: `r1.fq` and `r2.fq`

### B

Data folder: `/project/varcall_training/data/partial/fastq/pb/`
Files: `lr.fq`

## Suggested computational resources

To process the small test dataset, we suggest the following computational resources:

- CPUs: 1
- Memory: 8 GB

## Commands

### A

Navigate to your desired output folder and run the following command:

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	minimap2 \
	-a \
	-x sr \
	-t 1 \
	-R "@RG\tID:HG002\tPL:ILLUMINA\tSM:HG002_ILLUMINA" \
        /project/varcall_training/data/partial/genome/hg38/chr20.fa \
	/project/varcall_training/data/partial/fastq/illumina/r1.fq \
	/project/varcall_training/data/partial/fastq/illumina/r2.fq > 
	sr.sam 
```

Generate .bam, sort and index

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	samtools \
        sort \
	--write-index \
        sr.sam > sr.bam
```

### B

Navigate to your desired output folder and run the following command:

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	minimap2 \
        -a \
        -x map-hifi \
	-t 1 \
	-R "@RG\tID:HG002\tPL:PACBIO\tSM:HG002_PACBIO" \
	/project/varcall_training/data/partial/genome/hg19/20.fa \
        /project/varcall_training/data/partial/fastq/pb/lr.fq > lr.sam 
```

Generate .bam, sort and	index

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
        /project/varcall_training/bin/varcall_latest.sif \
	samtools \
	sort \
        --write-index \
        lr.sam > lr.bam
```

### Options explained

## minimap2

| Option | Description |
|--------|-------------|	
| `-a` | Output .sam (output .paf otherwise) |
| `-x` | Presets (minimap2 --help for all presets) |
| `-t` | Threads |
| `-R` | Add read group information in .sam |
| `<in.ref>` | Reference .fasta |
| `<in.read1>` | Read 1 - (single-end or read 1 in paired end) |
| `<in.read2>` | Read 2 (in paired end mode) |

## samtools

| Option | Description |
|--------|-------------|
| `--write-index` | Write index on the fly |


### Some statistics

Check statistics for the generated .bam files with `samtools stats` or `samtools flagstats`!
