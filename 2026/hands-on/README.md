# Hands-on documentation

- [Hands-on documentation](#hands-on-documentation)
	- [Introduction](#introduction)
	- [Materials](#materials)
		- [Data](#data)
		- [Software](#software)
	- [How to run the commands](#how-to-run-the-commands)
		- [Suggested computational resources](#suggested-computational-resources)
	- [Workflow](#workflow)
	- [Instructors](#instructors)

## Introduction

This repository contains the hands-on documentation for the genomic variant calling training course at Human Technopole, coordinated by members of the [Genome Analysis Unit](https://github.com/HTGenomeAnalysisUnit). 

Starting from the provided FASTQ files, you can follow all the steps to perform the variant calling analysis. The documentation is organized in chapters, each one focusing on a specific step of the analysis.
Output of one step can be used in the next one, but we also provide pre-processed files to start from any step. See the documentation of each step for more details.

In general, you should execute all commands in a compute node in a dedicated folder you created inside your group space or in your home directory (see the [How to run the commands section](#how-to-run-the-commands) for more details).

## Materials

Please refer to the course [GitHub repository](https://github.com/HTGenomeAnalysisUnit/varcall_training/tree/main/2026) for the materials used in the training and instructions.

### Data 

The test datasets used in the training are available in the `/project/varcall_training/data/partial` folder on the HT cluster. The data is organized in subfolders according to the analysis step, the type of data and the technology used for sequencing.

### Software

All the tools used in the training are available in this Singularity container: `/project/varcall_training/bin/varcall_latest.sif` folder on the HT cluster.

## How to run the commands

Remember to execute all the command in a compute node. First of all, request an interactive session:

```bash
srun --pty -p cpu-interactive --mem 8G --cpus-per-task 1 --time 4:00:00 /bin/bash
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

## Instructors

- [Edoardo Giacopuzzi](https://humantechnopole.it/en/people/edoardo-giacopuzzi/)
- [Davide Bolognini](https://humantechnopole.it/en/people/davide-bolognini/)
- [Bruno Ariano](https://humantechnopole.it/en/people/bruno-ariano/)
