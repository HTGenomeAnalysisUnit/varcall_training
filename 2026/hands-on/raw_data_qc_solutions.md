# Raw Data QC: The solutions

## 🔬 Sample 1: The Baseline (The "Good" Sample)

Before hunting for errors, you must know what "good" looks like.

In the FastQC/fastp report:

- Check the "Per base sequence quality". (Should be Green/High)
- Check "Adapter Content". (Should be empty)
- Check "Sequence Length Distribution". (Should be a sharp peak at 150bp)

**✅ Action:** These reads are ready for alignment. No cleaning needed.

---

## 🔬 Sample 2: The adapters

Look specifically at the **Adapter Content** module in FastQC and/or fastp reports.

- *Question:* Which adapter sequence is detected?
- *Question:* At what cycle (roughly) does the adapter contamination start rising?

**The Fix:** These reads are shorter than the sequencer cycle (150bp), so the machine read into the plastic adapter. You must trim this off. Luckly `fastp` does this automatically for you.

---

## 🔬 Sample 3: The "Unmappable"

- FastQC might look surprisingly normal (or show a little weird GC distribution).
- Check the mapping statistics (`samtools flagstat`).
- What percentage of reads mapped to the human genome? For a good sample we expected >= 95% usually.

**The Fix:** This is likely **species contamination** (e.g., bacteria).

**Further action:** You may want to use a tool like `fastq_screen` or `kraken2` to identify the species, then filter the BAM file to keep only mapped reads.

---

## 🔬 Sample 4: The duplicates

- Check the **Sequence Duplication Levels** in FastQC report.
- Check the number and computed the proportion of duplicate reads in the `samtools flagstat` output.

You will see a high percentage of non-unique reads, especially compared to sample 1.

**Cause:** PCR over-amplification. It may be that the library complexity was low or too many PCR cycles were used. You may want to check with the wet-lab team.

**The Fix:** You generally do not "fix" the FASTQ. You align it first, then mark duplicates in the BAM file so the variant caller ignores them.

---

## 🔬 Sample 5: The Quality Drop

Look at the **Per base sequence quality** boxplots in FASTQC / fastp reports.

- The quality scores likely plummet (turn red/orange) toward the end of the read (3' end).

**Cause:** Reagent depletion or "phasing" on the sequencer.

**The Fix:** Apply **Sliding Window** trimming to cut off the bad tails. Luckly, `fastp` can do this for you.

---

## 🔬 Sample 6: The "N" Spike

Look at the **Per base N content** module in the FastQC / fastp reports.

- You will see a spike in "N" (uncalled bases) at specific cycles or throughout reads.

**Cause:** Air bubbles in the flowcell or camera focus errors.

**The Fix:** Filter out reads that have too many Ns. Again, `fastp` can do this for you.

**Command Hint:**

```bash
# fastp filters Ns by default (usually max 5 Ns allowed)
fastp -i input_R1.fq.gz -I input_R2.fq.gz ...

```

---

## 🔬 Sample 7: The Hidden Enemy (Advanced)

FastQC will likely show **Green/Pass** for everything. The data looks clean.

- You need to inspect the **Variant Allele Frequency (VAF)** or use a tool like `VerifyBamID`.
- In a pure sample, variants are 50% (Het) or 100% (Hom). In this sample, you will see a cluster of variants at ~10% frequency (the contamination).
- Contaminating reads likely lead to increased heterozygosity and unexpected VAF distributions.
- `VerifyBamID` results will report a high FREEMIX value (e.g., >3-5%).

**Cause:** This is Human-on-Human contamination (Sample Swapping).

**The Fix:** If contamination is high (>3-5%), the sample is usually discarded. If low, tools can sometimes adjust quality scores to account for it.

---

## 🔬 Sample 8: The Weird Lengths

Check **Sequence Length Distribution** in FastQC / fastp reports.

- Instead of a single peak at 150bp, you see a curve or plateau ranging from 50bp to 150bp.
- This data has **already been trimmed** (or pre-processed) before you got it.
- If the fraction of very short reads (<50bp) is high, or you see a peak at unexpected lengths, this may indicate problems during sequencing or library prep.

**The Fix:** Usually, no fix is needed. However, you must ensure your aligner (like BWA-MEM) handles variable read lengths correctly (most modern aligners do).

---

## 🛠 Command Cheat Sheet

| Task | Tool | Typical Command Structure |
| --- | --- | --- |
| **View QC** | FastQC | `fastqc sample_R1.fastq.gz` |
| **Trim Adapters/Qual** | fastp | `fastp -i R1.fq -I R2.fq -o outR1.fq -O outR2.fq --html report.html` |
| **Sample mix check** | VerifyBamID | `VerifyBamID --BamFile aln.bam --Reference genome.fa --SVDPrefix ref/variants.vcf` |
| **Check Mapping %** | Samtools | `samtools flagstat aln.bam > aln.flag-stats.txt` |
| **Check duplication rate** | Samtools | `samtools flagstat aln.bam > aln.flag-stats.txt` |
| **Inspect VAF** | bcftools | `bcftools stats variants.vcf > stats.txt` |
