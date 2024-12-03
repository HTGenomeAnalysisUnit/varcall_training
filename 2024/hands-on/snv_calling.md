# Variant calling

In this chapter we start from and aligned and recalibrated BAM file and we perform variant calling with gatk HaplotypeCaller [gatk tool](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller).

## Input data

Data folder: `/project/varcall_training/data/partial/alignment`

Files: `hg38.snps.bqsr.bam` and `hg38.snps.bqsr.bam.bai`

Data folder: `/project/varcall_training/data/partial/truth`

Files: `snp.ranges.grch38.bed`

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
	gatk --java-options "-Xmx7g" HaplotypeCaller \
	-I $PWD/hg38.snps.bqsr.bam \
    -R /project/varcall_training/data/partial/genome/hg38/chr20.fa \
    -ERC GVCF \
    -L /project/varcall_training/data/partial/truth/snp.ranges.grch38.bed \
    -O $PWD/hg38.snps.bqsr.g.vcf.gz
```

This will run for approximately 10 minutes using a single CPU and will create a single file in the current directory:

- `hg38.snps.bqsr.g.vcf.gz`: the gvcf containing the variant called only for a the sample HG002

### Step 2

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	gatk --java-options "-Xmx7g" CombineGVCFs \
	-V $PWD/hg38.snps.bqsr.g.vcf.gz \
    -V $PWD/HG003.snps.bqsr.g.vcf.gz \
    -R /project/varcall_training/data/partial/genome/hg38/chr20.fa \
    -L chr20  \
    -O $PWD/cohort.g.vcf.gz
```
This will run for few seconds only:

- `cohort.g.vcf.gz`: A VCF containing the pooled variant called from HG002 and HG003

### Step 3
```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	gatk --java-options "-Xmx7g" GenotypeGVCFs \
	-V $PWD/cohort.g.vcf.gz \
    -R /project/varcall_training/data/partial/genome/hg38/chr20.fa \
    -L chr20  \
    -O $PWD/snps.gatk.vcf.gz
```

This will run for few seconds only and will create a combined vcf file:

- `snps.gatk.vcf.gz`: The joint VCF file

### Step 4
```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	bcftools norm -m-any \
    --check-ref w \
    -f /project/varcall_training/data/partial/genome/hg38/chr20.fa \
	-O z -o $PWD/snps.norm.gatk.joint.vcf.gz $PWD/snps.gatk.vcf.gz

```

This will run for few seconds only and will create a combined vcf file:

- `snps.norm.gatk.joint.vcf.gz`: The joint VCF spliited by alleles

### Step5
```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	tabix -p vcf $PWD/snps.norm.gatk.joint.vcf.gz
```

This will run for few seconds only and will create an index for the vcf file:

- `snps.norm.gatk.joint.vcf.gz.tbi`: The index for the VCF


### Options explained

| Option | Description |
|--------|-------------|	
| `--java-options` | How many GB of memory we want to dedicate to this step |
| `-I` | Input bam file (need to be indexed)|
| `-R` | Reference genome |
| `-ERC` | The mode to create the GVCF file |
| `-L` | a bed file containing the position to map onto |
| `---known-sites` | set of known site used to build the recalibrator model |
| `-O` | Output file for the recal table or the bam file depending on which step you are running |
| `--bqsr-recal-file` | path for the recalibration table |

