# GATK Variant Calling: Genome in a Bottle Ashkenazi Trio Example

This README demonstrates how to perform variant calling using GATK (Genome Analysis Toolkit) with data from the Genome in a Bottle (GIAB) Ashkenazi Trio. The pipeline includes explanations of each step and common supplementary commands with `bcftools` and `samtools`.

---

## 1. About Genome in a Bottle and the Ashkenazi Trio

**Genome in a Bottle (GIAB)** is a project by NIST (National Institute of Standards and Technology) to provide reference genome materials and datasets for benchmarking genomic analyses. These resources help ensure accuracy, reproducibility, and transparency in genome sequencing pipelines.

**The Ashkenazi Trio** refers to a gold-standard GIAB dataset consisting of high-quality genome data from an Ashkenazi Jewish family: son (HG002), father (HG003), and mother (HG004). These benchmark samples are widely used to evaluate variant calling methods and pipelines.

### 1.1 Setting things up

Log to github codespaces and ask for the smallest machine with 2 cores.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/HTGenomeAnalysisUnit/varcall_training)

When you open the virtual machine you will see a VS code interface with terminal and you are located at `/workspaces/varcall_training`.
git clone https://github.com/HTGenomeAnalysisUnit/varcall_training.git

Once inside the repo follow these steps

cd varcall_training/2026/hands-on/varcall_day3/references

tar -xzvf chr20.fa.gz

cd ../varcall_day3

---

## 2. GATK Variant Calling Workflow

The following GATK commands process BAM files from the Ashkenazi Trio (HG002, HG003, HG004) to produce joint variant calls.

### 2.1. HaplotypeCaller in GVCF Mode

The following three commands run GATK HaplotypeCaller separately for each sample. This produces a GVCF file for each, capturing potential variant information per sample in a format suitable for joint genotyping.

```bash
mkdir gvcf_files

gatk --java-options "-Xmx7g" HaplotypeCaller \
  -I bam_files/HG004.chr20.subset.rg.smaller.bam \
  -R ../references/chr20.fa \
  -ERC GVCF \
  -L chr20_teaching_genes_10kbflank.bed \
  -O gvcf_files/HG004.g.vcf.gz 2> gvcf_files/HG004_varcall.log


gatk --java-options "-Xmx7g" HaplotypeCaller \
  -I bam_files/HG003.chr20.subset.rg.smaller.bam \
  -R ../references/chr20.fa \
  -ERC GVCF \
  -L chr20_teaching_genes_10kbflank.bed \
  -O gvcf_files/HG003.g.vcf.gz 2> gvcf_files/HG003_varcall.log

gatk --java-options "-Xmx7g" HaplotypeCaller \
  -I bam_files/HG002.chr20.subset.rg.smaller.bam \
  -R ../references/chr20.fa \
  -ERC GVCF \
  -L chr20_teaching_genes_10kbflank.bed \
  -O gvcf_files/HG002.g.vcf.gz 2> gvcf_files/HG002_varcall.log
```
_Do this also for `HG002` and `HG003` respectively by changing the input BAM and output GVCF filenames._

#### Command parameters explained:
- `--java-options "-Xmx7g"`: Allocate 7GB RAM to the Java process.
- `HaplotypeCaller`: GATK tool to call variants per sample.
- `-I <in.bam>`: Input BAM file with alignments (`HG004.chr20.subset.rg.bam`, etc).
- `-R <reference.fa>`: Reference genome fasta (`chr20.fa` covers just chromosome 20).
- `-ERC GVCF`: Emit reference confidence scores in GVCF format (required for joint genotyping).
- `-L <regions.bed>`: Limit calling to regions in the BED file (`chr20_teaching_genes.bed`).
- `-O <out.g.vcf.gz>`: Output GVCF file name.
- `2> <logfile>`: Write log output (stderr) to file.


#### View logs produced

### 2.2. CombineGVCFs

After producing GVCFs for each sample, combine them for joint genotyping.

```bash
mkdir combined_vcfs

gatk --java-options "-Xmx7g" CombineGVCFs \
  -V gvcf_files/HG002.g.vcf.gz \
  -V gvcf_files/HG003.g.vcf.gz \
  -V gvcf_files/HG004.g.vcf.gz \
  -R ../references/chr20.fa \
  -O combined_vcfs/cohort.vcf.gz 2> combined_vcfs/combine.log
```

**Purpose:**  
Combines multiple GVCFs into a single multi-sample GVCF.  

**Parameters:**
- `-V <gvcf>`: (Repeat for each sample) Input GVCF files.
- `-R <reference.fa>`: Reference genome file.
- `-O <out.vcf.gz>`: Output combined GVCF file.
- `2> combine.log`: Log file.

### 2.3. GenotypeGVCFs

Genotype the combined GVCF to produce classic, fully-genotyped VCF.

```bash
gatk --java-options "-Xmx7g" GenotypeGVCFs \
  -R ../references/chr20.fa \
  -V combined_vcfs/cohort.vcf.gz \
  -O combined_vcfs/cohort.joint.vcf.gz
```

**Purpose**:  
Takes the combined multi-sample GVCF, performs joint genotyping, and outputs a standard VCF listing variants and genotypes for each individual.

**Parameters:**
- `-R <reference.fa>`: Reference genome.
- `-V <cohort.vcf.gz>`: Combined GVCF input.
- `-O <out.vcf.gz>`: Output genotyped VCF.

---

## 3. About `bcftools`

**bcftools** is a versatile set of command-line utilities for manipulating and analyzing VCF (Variant Call Format) and BCF (Binary VCF) files, which are the standard file formats for storing genomic variant data.

**Key features:**
- Efficient viewing and filtering of variants.
- Powerful querying of specific fields or regions.
- Format conversion (VCF ⇔ BCF).
- Annotation, merging, and statistical analysis.
- Supports indexed access for fast region queries.

### Common Applications

1. **Viewing and Inspecting Variants:**
   ```bash
   bcftools view combined_vcfs/cohort.joint.vcf.gz | less

   bcftools view -h combined_vcfs/cohort.joint.vcf.gz | less

   bcftools view -H combined_vcfs/cohort.joint.vcf.gz | less

   bcftools stats combined_vcfs/cohort.joint.vcf.gz | less
   ```

   Check the stats (INFO across samples, FORMAT per sample)

2. **Filtering Data:**
   - By region (chr20:10000-20000):
     ```bash
     bcftools view -r chr20:50000000-60000000 combined_vcfs/cohort.joint.vcf.gz
     ```
   - By minimum quality (QUAL ≥30):
     ```bash
     bcftools view -i 'QUAL>=30' combined_vcfs/cohort.joint.vcf.gz

     bcftools view -i 'FORMAT/DP[0]>5 && FORMAT/DP[1]>0' combined_vcfs/cohort.joint.vcf.gz

     bcftools view -i 'GT[0]=="0/0"' combined_vcfs/cohort.joint.vcf.gz
     ```
   - Output results:
   ```bash
     bcftools view -i 'FORMAT/DP[0]>5 && FORMAT/DP[1]>0' combined_vcfs/cohort.joint.vcf.gz -Oz -o test_output_combined.vcf.gz
   ```


3. **Querying Specific Fields:**
   - Extract basic info (chrom, pos, ref, alt):
     ```bash
     bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' cohort.joint.vcf.gz
     ```
   - Extract genotype info for samples:
     ```bash
     bcftools query -f '%CHROM\t%POS\t[%SAMPLE=%GT\t]\n' cohort.joint.vcf.gz
     ```

4. **Statistics and Summary:**
   ```bash
   bcftools stats cohort.joint.vcf.gz
   ```

### Filtering Syntax in `bcftools`

- `-i` (`--include`) : Keep records matching an expression.
- `-e` (`--exclude`): Remove records matching an expression.
- Expression examples:
  - Quality: `-i 'QUAL>=30'`
  - Heterozygous variants: `-i 'GT="het"'`
  - SNPs only: `-i 'TYPE="snp"'`
  - Sample DP > 20: `-i 'FORMAT/DP>20'`
- Combine expressions with `&` (and), `|` (or), and parentheses.
Once you have your VCF files, `bcftools` is useful to inspect, filter, and analyze variants.

---

## 4. Using `samtools mpileup`: When and Why

`samtools mpileup` creates a pileup summary of read alignments per genomic position.

**Use cases include:**
- Visual inspection of read support at variant sites,
- Troubleshooting variant calls,
- Creating custom formats for downstream analyses.

**Example usage:**

```bash
samtools mpileup -f ../references/chr20.fa HG004.chr20.subset.rg.smaller.bam | less
```

### Parameters:
- `-f <reference.fa>`: Reference genome Fasta.
- `<file.bam>`: Input BAM alignment.
- `| less`: View pileup in paginated format.

_You might use mpileup to investigate whether reads really support the genotype at a specific site or to manually inspect difficult variants._

---


## Interpreting the `samtools mpileup` BASES Column

When examining the output from `samtools mpileup`, the "BASES" column encodes how bases from aligned reads relate to the reference sequence at each position. Here are some common symbols and what they mean:

| Symbol       | Meaning                                               |
|--------------|------------------------------------------------------|
| `.`          | Base matches the reference on the forward strand      |
| `,`          | Base matches the reference on the reverse strand      |
| `A C G T`    | Mismatch (ALT base) on the forward strand            |
| `a c g t`    | Mismatch (ALT base) on the reverse strand            |
| `^`          | Start of a read                                      |
| `^g`         | Start of read; next character (`g`) encodes mapping quality |
| `$`          | End of a read                                        |

### Example

Given the following pileup sequence for the BASES column:

```
..^g.
```

**Means:**

- `.` : reference match
- `.` : reference match
- `^g.` :
  - `^g`: a read starts at this position, where `g` represents the encoded mapping quality (ASCII-encoded, not always literally the letter 'g')
  - `.` : the base at this position matches the reference

So in this example, there are two reads that match the reference, then a new read starts (`^g`), and the first base of this new read also matches the reference.


## References

- [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)
- [GATK Documentation](https://gatk.broadinstitute.org/)
- [bcftools Documentation](http://samtools.github.io/bcftools/)
- [samtools mpileup](http://www.htslib.org/doc/samtools-mpileup.html)
