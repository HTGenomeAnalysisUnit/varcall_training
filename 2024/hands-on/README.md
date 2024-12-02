# Hands-on documentation

- [Hands-on documentation](#hands-on-documentation)
	- [Introduction](#introduction)
	- [Materials](#materials)
		- [Data](#data)
		- [Software](#software)
	- [How to run the commands](#how-to-run-the-commands)
	- [1. FASTQ QC](#1-fastq-qc)
	- [X. Variant annotation](#x-variant-annotation)

## Introduction

## Materials

Please refer to the course [GitHub repository](https://github.com/HTGenomeAnalysisUnit/varcall_training/tree/main/2024) for the materials used in the training and instructions.

### Data 

The test datasets used in the training are available in the `/project/varcall_training/data/partial` folder on the HT cluster. The data is organized in subfolders according to the analysis step, the type of data and the technology used for sequencing.

### Software

All the tools used in the training are available in this Singularity container: `/project/varcall_training/bin/varcall_latest.sif` folder on the HT cluster.

## How to run the commands

In general, to run the commands in the Singularity container, you can use the following syntax:

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	<command> <options>
```

## 1. FASTQ QC

The first step in our analysis is to perform quality control (QC) on the raw FASTQ files. We will use the `fastp` tool to perform this task. This task is performed only for short-reads data.

Follow the steps in the [FASTQ QC](fastq_qc.md) tutorial to perform the QC.

## X. Variant annotation

The last step in our analysis is to annotate the variants called in the previous steps. We will use:

- `bcftools csq` tool to annotate the VCF file with the consequence of the variants
- `vcfanno` tool to annotate the VCF file with additional information

Follow the steps in the [Variant annotation](variant_annotation.md) tutorial to perform the annotation.