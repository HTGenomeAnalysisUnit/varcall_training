# Raw data Quality Control & Forensics Lab

In this chapter we start from raw FASTQ files of a short read sequencing experiment and we perform quality control (QC) and trimming of the reads using [fastp tool](https://github.com/OpenGene/fastp).

## Suggested computational resources

To process the small test dataset, we suggest the following computational resources:

- CPUs: 2
- Memory: 12.GB

Remember to execute all the command in a compute node. First of all, request an interactive session:

```bash
srun --pty -p cpu-interactive --mem 12G --cpus-per-task 2 --time 4:00:00 /bin/bash
```

## Objective

You have received sequencing data for Human Chromosome 20 from 8 different sequencing runs performed on Illumina sequencer using a 2x150bp paired-end protocol.

Based on that you expect:

- High-quality reads (Phred score > 30 across most of the read length).
- Properly paired reads with mean size around 150bp.
- Low levels of uncalled bases (N content < 1%).
- High mapping rate to the human reference genome (>98%).
- Low level of sequence duplication (5-7% or less).

Your job is to use Quality Control (QC) tools to diagnose each dataset and identify any specific problem they may have, and propose a fix.

We will see in the next day how to align the reads and call variants, but for now we will focus on diagnosing and fixing the raw data quality issues.

---

## Tools used

- `fastp`: All-in-one pre-processor (adapters, quality trimming, N-filtering).
- `fastqc`: Inspect main reads quality metrics.
- `VerifyBamID`: Checks for sample contamination and swaps.
- `samtools`: Collect statistics from BAM files.
- `bcftools`: Compute variant statistics.

---

## 📂 The Datasets

You will find the data in the following directory:

```text
/project/varcall_training/data/simulated/exercise_data
```

Here you can find 8 different samples (`sample_1` to `sample_8`), each corresponding to a different sequencing scenario with specific issues to diagnose and fix. In each folder you can find a pair of FASTQ files and a BAM file.

You can also find a `merged_variants.vcf.gz` file containing the expected variants for all samples.

If you have trouble running the commands, the expected outputs are in the `/project/varcall_training/data/simulated/solution_data` folder.

---

## 🛠️ Step-by-Step Protocol

### 1. QC & Trimming (`fastp` and `fastqc`)

Before aligning, we must clean the data. `fastp` automatically detects adapters, trims low-quality bases (Phred < 15), and removes reads with too many Ns.

```bash
# Variables
SAMPLE_ID="sample1" # Change for each sample
R1="/project/varcall_training/data/simulated/exercise_data/${SAMPLE_ID}/${SAMPLE_ID}_R1.fastq.gz"
R2="/project/varcall_training/data/simulated/exercise_data/${SAMPLE_ID}/${SAMPLE_ID}_R2.fastq.gz"

# Run fastp
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
fastp \
	--thread 2 \
	--detect_adapter_for_pe \
	-i $R1 \
	-I $R2 \
	-o ${SAMPLE_ID}_clean_R1.fq.gz \
	-O ${SAMPLE_ID}_clean_R2.fq.gz \
	--json ${SAMPLE_ID}_fastp.json \
	--html ${SAMPLE_ID}_fastp.html

# Run fastqc
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
fastqc \
	--threads 2 \
	--outdir fastqc_${SAMPLE_ID} \
	$R1 \
	$R2
```

---

### 2. Mapping Statistics (`samtools flagstat`)

Now we generate a summary report to see how well our data aligned.

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
samtools flagstat /project/varcall_training/data/simulated/exercise_data/${SAMPLE_ID}/${SAMPLE_ID}.dedup.sort-coordinate.bam > ${SAMPLE_ID}_flagstat.txt
```

---

### 3. Compute sample mixing contamination (`VerifyBamID`)

We use `VerifyBamID` to estimate contamination levels in our sample.

Reference files to be used with VerifyBamID are already available in `/project/varcall_training/data/simulated/ref`.

In the real-life you will probably use files from the 1000 Genomes Project or gnomAD.

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
VerifyBamID \
  --SVDPrefix /project/varcall_training/data/simulated/ref/population_with_genotypes.vcf \
  --BamFile /project/varcall_training/data/simulated/exercise_data/${SAMPLE_ID}/${SAMPLE_ID}.dedup.sort-coordinate.bam \
  --Reference /project/varcall_training/data/simulated/ref/chr20.fa \
  --NumThread 4 \
  --Output ${SAMPLE_ID}_verifybamid
```

---

### 4. Compute variant statistics (`bcftools stats`)

You have called variants for the 8 samples in the file `merged_variants.vcf.gz` (see [dataset](#-the-datasets)). We will see how to do call variants later in the course, for now, we want to compute variant statistics using `bcftools stats`.

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
bcftools stats -s- /project/varcall_training/data/simulated/exercise_data/merged_variants.vcf.gz > merged_variants.stats
```

We can then use `plot-vcfstats` to visualize the statistics.

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/varcall_latest.sif \
plot-vcfstats -p vcf_plots merged_variants.stats
```

In the `vcf_plots` folder you will find multiple plots summarizing the variant statistics and a `summary.pdf` file with all plots combined. Expecially, you should look at 

- `snps_by_sample.0.png`: Number of variants per sample.
- `het_by_sample.0.png`: Heterozigosity. Number of heterozygous calls / number of homozigous calls per sample.

Keep in mind here the samples are numbered from 0 to 7 (not from 1 to 8).

---

## Inspect Logs & Reports

Now that you have successfully generated all the QC logs and reports for all sample, review the metrics to identify which issue may affect each sample.

### The `fastqc` and `fastp` HTML report

These reports contain multiple modules with different quality metrics, describing various aspects of the sequencing data quality. Look in particular at:

- **Insert Size:** Is the peak at 150-200bp? Do you see any abnormal peak or left tail (short inserts)?
- **Per Base Sequence Quality:** Are there low-quality tails at the start/end of reads?
- **Per Base Sequence Content:** Are there abnormal patterns (e.g., high uneven distribution of nucleotides or high N content)?
- **Adapter Trimmed:** Was there significant adapter contamination on the read end?

### The `samtools flagstat` text file

Run `cat ${SAMPLE_ID}_flagstat.txt` and check:

1. **Mapped %:** Should be >98% for a good human sample.
2. **Properly Paired %:** Should be high (>95%).
3. **Duplicates:** Look at the number of duplicated reads. A good WGS sample should have about 5-7% duplicated reads among mapped reads. (Note: `samblaster` marks them but doesn't remove them).

### The `VerifyBamID` text file

Look for the `FREEMIX` value in the `${SAMPLE_ID}_verifybamid.selfSM` file. This value estimates the fraction of contamination in the sample measuring the proportion of reads coming from a different individual. A good sample should have `FREEMIX < 0.03` (i.e., <3% contamination).

### The `bcftools stats` plots

You have to look at the PDF file `vcf_plots/summary.pdf` and in particular at the plots about number of variants per sample, number of het calls and ratio of het/hom.

Combine these observations with the VerifyBamID results and samtools flagstat to identify sample issues.
