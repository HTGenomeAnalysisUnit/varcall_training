# Hands-on documentation

- [Hands-on documentation](#hands-on-documentation)
	- [Introduction](#introduction)
	- [Materials](#materials)
		- [Data](#data)
		- [Software](#software)
	- [How to run the commands](#how-to-run-the-commands)
		- [Suggested computational resources](#suggested-computational-resources)
	- [Workflow](#workflow)
		- [1. FASTQ QC](#1-fastq-qc)
		- [2. Alignment](#2-alignment)
		- [3. Mark duplicates](#3-mark-duplicates)
		- [4. Base recalibration](#4-base-recalibration)
		- [5. Variant calling](#5-variant-calling)
		- [6. Variant annotation](#6-variant-annotation)
		- [7. Variant visualization](#7-variant-visualization)
		- [8. Variant benchmarking](#8-variant-benchmarking)
	- [Instructors](#instructors)

## Introduction

This repository contains the hands-on documentation for the genomic variant calling training course at Human Technopole, coordinated by members of the [Genome Analysis Unit](https://github.com/HTGenomeAnalysisUnit). 

Starting from the provided FASTQ files, you can follow all the steps to perform the variant calling analysis. The documentation is organized in chapters, each one focusing on a specific step of the analysis.
Output of one step can be used in the next one, but we also provide pre-processed files to start from any step. See the documentation of each step for more details.

In general, you should execute all commands in a compute node in a dedicated folder you created inside your group space or in your home directory (see the [How to run the commands section](#how-to-run-the-commands) for more details).

## Materials

Please refer to the course [GitHub repository](https://github.com/HTGenomeAnalysisUnit/varcall_training/tree/main/2024) for the materials used in the training and instructions.

### Data 

The test datasets used in the training are available in the `/project/varcall_training/data/partial` folder on the HT cluster. The data is organized in subfolders according to the analysis step, the type of data and the technology used for sequencing.

### Software

All the tools used in the training are available in this Singularity container: `/project/varcall_training/bin/varcall_latest.sif` folder on the HT cluster.

## How to run the commands

Remember to execute all the command in a compute node. First of all, request an interactive session:

```bash
srun --pty -p interactive --mem 8G --cpus-per-task 1 --time 4:00:00 /bin/bash
```

You should created a dedicated folder for the training in your group space or in your home directory. For example you can create a varcall_training folder in you home directory

```bash
mkdir ~/varcall_training
cd ~/varcall_training
```

Then you need to load the singularity module:

```bash
module load singularity
```

In general, to run the commands in the Singularity container, you can use the following syntax:

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	<command> <options>
```

The `-B` option is used to bind the directories to the container. The `$PWD` variable is used to bind the current working directory to the container. The `/project/varcall_training` directory is used to bind the data folder containig all the training material. The `/localscratch` directory is used to bind the scratch directory of the compute node which contains the temp folder.

### Suggested computational resources

To process the small test dataset, we suggest the following computational resources for most of the steps:

- CPUs: 1
- Memory: 8 GB

## Workflow

### 1. FASTQ QC

The first step in our analysis is to perform quality control (QC) on the raw FASTQ files. We will use the `fastp` tool to perform this task. This task is performed only for short-reads data.

Follow the steps in the [FASTQ QC](fastq_qc.md) tutorial to perform the QC.

### 2. Alignment

The second step in our analysis is to align the reads to the reference genome and perform some QC. We will use the `minimap2` tool to perform this task.

Follow the steps in the [Alignment](alignment.md) tutorial to perform the alignment.

### 3. Mark duplicates

The third step in our analysis is to mark the duplicates in the aligned reads. We will use the `samblaster` tool to perform this task. This task is performed only for short-reads data.

Follow the steps in the [Mark duplicates](mark_duplicates.md) tutorial to mark the duplicates.

### 4. Base recalibration

### 5. Variant calling

We will now perform small and structural variant calling.

Follow the steps in the [SV calling](sv_calling.md) tutorial to call SVs in your alignment.

### 6. Variant annotation

The last step in our analysis is to annotate the variants called in the previous steps. 

For small variants, we will use:

- `bcftools csq` tool to annotate the VCF file with the consequence of the variants
- `vcfanno` tool to annotate the VCF file with additional information

For structural variants, we will use:

- `SVAfotate` tool to annotate the VCF file with population allele frequencies from 1000G and gnomAD databases

Follow the steps in the [Variant annotation](variant_annotation.md) tutorial to perform the annotation.

### 7. Variant visualization

### 8. Variant benchmarking

(Optional). We can check how our variant callers performed when comparing to a trugh set.

Follow the steps in the [Variant benchmark](benchmark.md) tutorial to perform the benchmarking.

## Instructors

- [Edoardo Giacopuzzi](https://humantechnopole.it/en/people/edoardo-giacopuzzi/)
- [Davide Bolognini](https://humantechnopole.it/en/people/davide-bolognini/)
- [Bruno Ariano](https://humantechnopole.it/en/people/bruno-ariano/)

Seminars

- [Andrea Guarracino](https://andreaguarracino.github.io/)
- [Alessandro Raveane](https://humantechnopole.it/en/people/alessandro-raveane/)
- [Chela Tandiwe James](https://humantechnopole.it/en/people/chela-tandiwe-james/)
