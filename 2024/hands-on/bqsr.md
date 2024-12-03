# Base score quality recalibration (BSQR)

In this chapter we start from aligned BAM files of a short read sequencing experiment and we perform base call quality recalibration (BQSR) [gatk tool](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR).

## Input data

Data folder: `/project/varcall_training/data/partial/alignment`

Files: `hg38.snps.bam` and `hg38.snps.bam.bai`

## Suggested computational resources

To process the small test dataset, we suggest the following computational resources:

- CPUs: 1
- Memory: 8 GB

## Commands

Navigate to your desired output folder and run the following command:

### Step 1 

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	gatk --java-options "-Xmx7g" BaseRecalibrator \
	-I /project/varcall_training/data/partial/alignment/hg38.snps.bam \
    -R /project/varcall_training/data/partial/genome/hg38/chr20.fa \
    --known-sites /project/varcall_training/data/partial/annotations/snps/dbsnp_146.hg38.vcf.gz \
     -O $PWD/recal_data.table
```

This will run for approximately 3 minutes using a single CPU and will create a single file in the current directory:

- `recal_data.table`: the table used for step 2 for recalibrating the variants

### Step 2

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	gatk --java-options "-Xmx7g" ApplyBQSR \
	-I /project/varcall_training/data/partial/alignment/hg38.snps.bam \
    -R /project/varcall_training/data/partial/genome/hg38/chr20.fa \
    --bqsr-recal-file $PWD/recal_data.table \
    -O $PWD/hg38.snps.bqsr.bam
```
This will run for approximately 2 minutes using a single CPU and will create a single BAM file in the current directory:

- `hg38.snps.bqsr.bam`: the table used for step 2 for recalibrating the variants

### Step 3

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	samtools index $PWD/hg38.snps.bqsr.bam
```
This will run for few seconds only and will create an index for the BAM file:

- `hg38.snps.bqsr.bam.bai`

### Options explained

| Option | Description |
|--------|-------------|	
| `--java-options` | How many GB of memory we want to dedicate to this step |
| `-I` | Input bam file (need to be indexed)|
| `-R` | Reference genome |
| `---known-sites` | set of known site used to build the recalibrator model |
| `-O` | Output file for the recal table or the bam file depending on which step you are running |
| `--bqsr-recal-file` | path for the recalibration table |
