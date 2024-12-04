# Benchmark

We can check how our variant callers performed comparing variant calls to a ground trugh using [truvari](https://github.com/ACEnglish/truvari)

## Input data

### A

Truth folder: `/project/varcall_training/data/partial/truth/snps/`
Files: `snp.vcf.gz`, `snp.vcf.gz.tbi`
Variants folder: `/project/varcall_training/data/partial/vcf`
Files: `snps.norm.gatk.vcf.gz`, `snps.norm.gatk.vcf.gz.tbi`

### B

Truth folder: `/project/varcall_training/data/partial/fastq/truth/svs/`
Files: `svs.vcf.gz`, `svs.vcf.gz.tbi`
Variants folder: `/project/varcall_training/data/partial/vcf`
Files: `svs.vcf.gz`, `svs.vcf.gz.tbi`

## Suggested computational resources

To process the small test dataset, we suggest the following computational resources:

- CPUs: 1
- Memory: 1 GB

## Commands

### A

Navigate to your desired output folder and run the following command:

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	truvari \
	bench \
	-s 0 \
	-b /project/varcall_training/data/partial/truth/snps/snp.vcf.gz \
	-c /project/varcall_training/data/partial/vcf/snps.norm.gatk.vcf.gz \
	-o snps_benchmark
```

### B

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
        /project/varcall_training/bin/varcall_latest.sif \
        truvari	\
        bench \
	-b /project/varcall_training/data/partial/truth/svs/svs.vcf.gz \
        -c /project/varcall_training/data/partial/vcf/svs.vcf.gz \
        -o svs_benchmark
```

### Options explained

| Option | Description |
|--------|-------------|	
| `-b` | Ground truth vcf (base vcf) |
| `-c` | Vcf from variant caller (comparison vcf) |
| `-o` | Output folder |
| `-s` | Minimum variant size to consider from -c|

### Some statistics

Check a couple of lines from the generated vcf files in the output folder with `bcftools view`
