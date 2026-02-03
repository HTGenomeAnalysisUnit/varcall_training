# Variant annotation

In this chapter we start from QCed variants in VCF format and we perform annotation of the variants using `bcftools csq` and `vcfanno` tools for small variants and `SVAFotate` for structural variants.

- [Variant annotation](#variant-annotation)
	- [Suggested computational resources](#suggested-computational-resources)
	- [Small variants](#small-variants)
		- [Input data](#input-data)
		- [Commands](#commands)
			- [1. Annotate gene consequences](#1-annotate-gene-consequences)
			- [2. Annotate additional information](#2-annotate-additional-information)
	- [Structural variants](#structural-variants)
		- [Input data](#input-data-1)
		- [Commands](#commands-1)

## Suggested computational resources

To process the small test dataset, we suggest the following computational resources:

- CPUs: 1
- Memory: 8 GB

## Small variants

For small variant we will first use [bcftools csq](https://samtools.github.io/bcftools/howtos/csq-calling.html) to annotate gene consequences and then we will use [vcfanno](https://github.com/brentp/vcfanno) to add a set of useful annotations to our input VCF.

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

This will run for approximately 1 minute using a single CPU and will create a new vcf file named `annotated_csq.vcf.gz` and the corresponding index file in the current directory. This file contains annotations of gene consequences.

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

For structural variants we will use [SVAFotate](https://github.com/fakedrtom/SVAFotate) to annotate each variant with allele frequencies in different populations and the fraction of each variant never observed in the external populations.

### Input data

Data folder: `/project/varcall_training/data/partial/vcf/`

Input VCF file: `svs.vcf.gz`

### Commands

```bash
singularity exec -B $PWD -B /localscratch -B /project/varcall_training \
	/project/varcall_training/bin/varcall_latest.sif \
	svafotate annotate --cpu 1 -O vcfgz -l 500000 -c 0.001 -f 0.5 \
	-b /project/varcall_training/data/partial/annotations/svs/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz \
	-v /project/varcall_training/data/partial/vcf/svs.vcf.gz \
	-o annotated_sv.vcf.gz
```

This will run for approximately 3 minutes using a single CPU and will create a new vcf file named `annotated_sv.vcf.gz` in the current directory. This file contains a set of new INFO fields representing allele frequencies of the structural variants in different populations (1000G, gnomAD, etc.) and the fraction of each variant never observed in the external populations.

**Options explained**

| Option | Description |
|--------|-------------|
| `--cpu` | Number of threads to use |
| `-O` | Output format |
| `-l` | Maximum length of the SVs in annotation source to be included. This will avoid to mess up annotation by including very large events that will overlap with everything |
| `-c` | Minimum allele frequency to consider a variant when computing the fraction of SVs never observed |
| `-f` | Fraction of reciprocal overlap to consider 2 events of the same type the same |
| `-b` | Path to the bed file with allele frequencies |
| `-v` | Path to the input VCF file |
| `-o` | Output file path |