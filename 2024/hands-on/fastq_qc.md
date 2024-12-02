# FASTQ qc

In this chapter we start from raw FASTQ files of a short read sequencing experiment and we perform quality control (QC) and trimming of the reads using [fastp tool](https://github.com/OpenGene/fastp).

## Input data

Data folder: `/project/varcall_training/data/partial/fastq/illumina`

Files: `r1.fq` and `r2.fq`

## Suggested computational resources

To process the small test dataset, we suggest the following computational resources:

- CPUs: 1
- Memory: 8 GB

## Commands

Navigate to your desired output folder and run the following command:

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
	fastp \
	-i /project/varcall_training/data/partial/fastq/illumina/r1.fq \
	-I /project/varcall_training/data/partial/fastq/illumina/r2.fq \
	-o $PWD/r1.fastp.fq.gz -O $PWD/r2.fastp.fq \
	-h fastp.html -j fastp.json \
	-l 30 -r -p -w 1 -V
```

This will run for approximately 2 minutes using a single CPU and will create a set of files in the current directory:

- `r1.fastp.fq.gz` and `r2.fastp.fq.gz`: the trimmed and filtered FASTQ files
- `fastp.html`: a report of the QC performed by fastp
- `fastp.json`: a JSON file with the QC metrics

You can open the HTML report in a browser to inspect the QC results. 

The trimmed and filtered FASTQ files can be used as input for alignement and variant calling.

### Options explained

| Option | Description |
|--------|-------------|	
| `-i` | Input file for read 1 |
| `-I` | Input file for read 2 |
| `-o` | Output file for read 1 |
| `-O` | Output file for read 2 |
| `-h` | Output file for HTML report |
| `-j` | Output file for JSON report |
| `-l` | Minimum read length after trimming. Reads shorter than these will be discarded |
| `-r` | Activate sliding window quality trimming.<br>Move a sliding window from front to tail, if meet one window with mean quality < threshold,<br>drop the bases in the window and the right part, and then stop. |
| `-p` | Enable overrepresented sequence analysis |
| `-w` | Number of threads to use (asjust to match request n CPUs) |
| `-V` | Verbose log |