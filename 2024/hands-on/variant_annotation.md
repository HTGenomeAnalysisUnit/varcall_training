# Variant annotation

In this chapter we start from QCed variants in VCF format and we perform annotation of the variants using `bcftools csq` and `vcfanno` tools for small variants and `SVAFotate` for structural variants.

## Suggested computational resources

To process the small test dataset, we suggest the following computational resources:

- CPUs: 1
- Memory: 8 GB

## Small variants

For small variant we will first use `bcftools csq` to annotate gene consequences and then we will use `vcfanno` to add a set of useful annotations to our input VCF.

### Input data

Data folder: `/project/varcall_training/data/partial/vcf/`

Input VCF file: `snps.norm.bcftools.vcf.gz`

### Commands

#### 1. Annotate gene consequences

Navigate to your desired output folder and run the following command:

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	bcftools csq -l --write-index \
	-g /project/varcall_training/data/partial/annotations/gene_consequences/Homo_sapiens.GRCh38.109.gff3.gz \
	-f /project/varcall_training/data/partial/genome/hg38/chr20.fa \
	-Oz -o annotated_csq.vcf.gz \
	/project/varcall_training/data/partial/vcf/snps.norm.bcftools.vcf.gz
```

This will run for approximately 1 minute using a single CPU and will create a new vcf file named `annotated_csq.vcf.gz` in the current directory. This file contains annotations of gene consequences.

**Options explained**

| Option | Description |
|--------|-------------|
| `-l`   | Disable phased consequences since we don't have phased variants |
| `-f`   | Path to reference genome fasta file |
| `-g`   | Path to a gff3 file describing transcript structure |
| `-Oz`  | Set output formate to bgzip compressed VCF |
| `-o`   | Output file name |

#### 2. Annotate additional information

The annotations sources and fields are defined in the toml configuration file: `/project/varcall_training/data/partial/annotations/snps/annotations.toml`

You can now use the file you just generated to add some more annotations using `vcfanno`

Navigate to your desired output folder and run the following command:

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	vcfanno -p 1 -base-path /project/varcall_training/data/partial/annotations/snps \
	/project/varcall_training/data/partial/annotations/snps/annotations.toml \
	annotated_csq.vcf.gz > annotated.vcf
```

This will run for approximately 1 minute using a single CPU and will create a new vcf file named `annotated.vcf` in the current directory. This file contains a set of new INFO fields representing the annotations described in the `annotations.toml` file.

**Options explained**

| Option | Description |
|--------|-------------|
| `-p` | Number of threads to use |
| `-base-path` | Base path for the annotation files. This path is prefixed to all files defined in the toml configuration |

## Structural variants

For s