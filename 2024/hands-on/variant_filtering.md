# Variant recalibration and filtering (VQSR)

In this chapter we start from a joint VCF called file and we perform variant quality recalibration and filtering (VQSR) [gatk tool](https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator).

## Input data

Data folder: `/project/varcall_training/data/partial/vcf`

Files: `snps.norm.gatk.joint.vcf.gz` and `snps.norm.gatk.joint.vcf.gz.tbi`

Data folder: `/project/varcall_training/data/partial/gatk_reference`

Files: `dbsnp_146.hg38.vcf.gz` and `hapmap_3.3.hg38.vcf.gz`

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
	gatk --java-options "-Xmx7g" VariantRecalibrator \
    -V /project/varcall_training/data/partial/vcf/snps.norm.gatk.joint.vcf.gz \
    --trust-all-polymorphic \
    -mode BOTH  \
    --max-gaussians 1 \
    --resource:hapmap,known=false,training=true,truth=true,prior=15 /project/varcall_training/data/partial/gatk_reference/hapmap_3.3.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=7 /project/varcall_training/data/partial/gatk_reference/dbsnp_146.hg38.vcf.gz \
    -an QD \
    -an MQRankSum \
    -an ReadPosRankSum \
    -an FS \
    -an MQ \
    -an SOR \
    -an DP \
    -O $PWD/cohort_snps.recal  \
    --tranches-file $PWD/cohort_snps.tranches
```


This will run for only few seconds using a single CPU and will create a single file in the current directory:

- `cohort_snps.recal`: the table used for step 2 for recalibrating the variants
- `cohort_snps.tranches`: the table with the tranches where each variant belong

### Step 2

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	gatk --java-options "-Xmx7g" ApplyVQSR \
	-V /project/varcall_training/data/partial/vcf/snps.norm.gatk.joint.vcf.gz \
    -R /project/varcall_training/data/partial/genome/hg38/chr20.fa \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file $PWD/cohort_snps.tranches \
    --recal-file $PWD/cohort_snps.recal \
    -O $PWD/snps.norm.gatk.joint.vsqr.vcf.gz
```
This will run for few seconds only using a single CPU and will create a single VCF file in the corrected variant score:

- `snps.norm.gatk.joint.vsqr.vcf.gz`: the vcf with the corrected variant quality score
- `snps.norm.gatk.joint.vsqr.vcf.gz.tbi`: the index for the vcf

### Options explained

| Option | Description |
|--------|-------------|	
| `--java-options` | How many GB of memory we want to dedicate to this step |
| `-V` | A variant calling file VCF (need to be indexed) |
| `-R` | Reference genome |
| `--trust-all-polymorphic` | Trust that all the input training sets' unfiltered records contain only polymorphic sites to drastically speed up the computation |
| `-O` | Output file for the recal table or the VCF file depending on which step you are running |
| `mode` | Either SNP, INDEL or BOTH, depending what you want to recalibrate |
| `resource` | Which reference to use to train and validate the ML model|
| `an` | Which annotation to use to train the model|
| `--truth-sensitivity-filter-level` | Where to set the sensitivity filter|
