# Hands-on documentation

- [Hands-on documentation](#hands-on-documentation)
	- [Introduction](#introduction)
	- [Materials](#materials)
		- [Data](#data)
		- [Software](#software)
	- [How to run the commands](#how-to-run-the-commands)
	- [1. FASTQ QC](#1-fastq-qc)
	- [2. Alignment](#2-alignment)
	- [3. Mark duplicates](#3-mark-duplicates)
	- [4. Base recalibration](#4-base-recalibration)
	- [5. Variant calling](#5-variant-calling)
	- [6. Variant annotation](#6-variant-annotation)
	- [7. Variant visualization](#7-variant-visualization)

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

## 2. Alignment

The second step in our analysis is to align the reads to the reference genome. We will use the `minimap2` tool to perform this task.

Follow the steps in the [Alignment](alignment.md) tutorial to perform the alignment.

## 3. Mark duplicates

The third step in our analysis is to mark the duplicates in the aligned reads. We will use the `samblaster` tool to perform this task. This task is performed only for short-reads data.

Follow the steps in the [Mark duplicates](mark_duplicates.md) tutorial to mark the duplicates.

## 4. Base recalibration

## 5. Variant calling

## 6. Variant annotation

The last step in our analysis is to annotate the variants called in the previous steps. 

For small variants, we will use:

- `bcftools csq` tool to annotate the VCF file with the consequence of the variants
- `vcfanno` tool to annotate the VCF file with additional information

For structural variants, we will use:

- `SVAfotate` tool to annotate the VCF file with population allele frequencies from 1000G and gnomAD databases

Follow the steps in the [Variant annotation](variant_annotation.md) tutorial to perform the annotation.

## 7. Variant visualization
