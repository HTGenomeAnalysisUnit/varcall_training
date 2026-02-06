## Codespace

As in the previous days, you can access all the material in the [codespace](https://github.com/codespaces/new?skip_quickstart=true&machine=basicLinux32gb&repo=883029421&ref=main&devcontainer_path=.devcontainer%2Fdevcontainer.json&geo=EuropeWest) for this repo.

# Pangenome practical session

We will build a pangenome graph for locus CYP2D6-CYP2D7 and use it to genotype a short read sample with COSIGT. Additionally, we will convert the graph into a vcf with vg toolkit.

## 1. Input files

Input file `haplotypes_plus_ref.fa` contains the haplotype sequences from 8 genomes plus the reference genome GRCh38. We already extracted the relevant region of each genome assembly using `minimap2` and `impg`.

```bash
samtools faidx haplotypes_plus_ref.fa
cut -f 1,2 haplotypes_plus_ref.fai
```

These haplotypes were taken from HGSVCv3 and the names follow the PanSN specification:

SampleName#haplNumber#contigName

## 2. Building the graph

We will use `pggb` to build the graph. Then we will extract the name of the reference path with `odgi paths` and sort the graph relative to the reference. 

```bash
pggb -c 2 -k 101 -x auto -i haplotypes_plus_ref.fa -n 9 -o cyp2_graph
odgi paths -i cyp2_graph/*smooth.final.og -L | grep grch > ref_path.txt  
odgi sort -Y -H ref_path.txt  -i cyp2_graph/*smooth.final.og -o cyp2_graph/cyp.og
```
## 3. Inspection

Within the folder `cyp2_graph` there are all of `pggb`'s output files. In particular there will be a png ending with `.viz_depth_multiqc.png`. This is a heatmap of node coverage across the entire graph. 

## 4. Parsing

Since the `.og` format is not human-readable we can use odgi view to convert it to `gfa`. We can parse this gfa file to extract node length with `awk` by counting the number of characters on "S" lines. We can also use `odgi paths` to get information on the paths within the graph. For example, `odgi paths -H` will print the node coverage of each haplotype in a tab-delimited format (`path.name`, `path.length`, `path.step.count`, `node.1`, `node.2`, ...). `odgi paths -L` which we used before will instead print the path names only.

```bash
odgi view -i cyp2_graph/cyp.og -g > cyp2.gfa
zgrep '^S' cyp2.gfa | awk '{{print("node."$2,length($3))}}' OFS="\\t" > cyp2.node_length.tsv
odgi paths -i cyp2.gfa -H | cut -f 1,4- | gzip > cyp2.node_cov.tsv.gz
```

## 5. Short Read Alignment and Genotyping

After building the pangenome graph, we can genotype a short-read sample against it using `cosigt`. This workflow involves aligning reads to the haplotype sequences, projecting those alignments onto the graph, and computing node coverage to identify the most similar haplotypes.

### 5.1 Indexing and Aligning Reads

First, we index the haplotype sequences with `bwa-mem2` and align the sample reads. The `-p` flag treats interleaved input, and `-h 1000` reports up to `1000` alternative hits to capture multi-mapping reads, which is crucial for accurate graph-based genotyping.

```bash
bwa-mem2 index haplotypes_plus_ref.fa
bwa-mem2 mem -t 2 -p -h 1000 haplotypes_plus_ref.fa na19983.fa.gz | samtools sort -T tmp -o na19983.realigned.bam --write-index
```

### 5.2 Projecting Alignments to the Graph

We use `gfainject` to convert the `bam` alignments to `gaf` format by mapping them to the reference paths in the `gfa` graph. The `--alt-hits 10000` parameter ensures we capture all alignment positions across the graph.
​
```bash
gfainject --gfa cyp2.gfa --bam na19983.realigned.bam --alt-hits 10000 | gzip > na19983.gaf
```

### 5.2 Computing Node Coverage

`gafpack` calculates node coverage from the GAF alignments. The `--len-scale` flag normalizes coverage by node length, and `--weight-queries` accounts for query occurrence frequencies to produce accurate coverage profiles.

```bash
gafpack --gfa cyp2.gfa --gaf na19983.gaf --len-scale --weight-queries | gzip > na19983.gafpack
```

### 5.3 Genotyping with COSIGT

`cosigt` uses cosine similarity between the sample's node coverage profile and the reference haplotype coverage to assign structural haplotypes. The `-p` parameter provides the reference node coverage, `-g` the sample coverage, and `-i` specifies the sample identifier.

```bash
cosigt -p cyp2.node_cov.tsv.gz -g na19983.gafpack -o na19983_genotype -i na19983
```

### 5.4 Masking Complex Regions

Complex or repetitive regions can reduce genotyping accuracy. We use `panplexity` to identify and mask such regions based on linguistic complexity metrics. The tool analyzes k-mer patterns across sliding windows to detect low-complexity sequences.

```bash
panplexity --input-gfa cyp2.gfa -t auto -k 16 -w 100 -d 100 --complexity linguistic -m cyp2.mask.txt --threads 2
```
We then re-run COSIGT with the mask file to exclude problematic regions from genotyping:

```bash
cosigt -p cyp2.node_cov.tsv.gz -g na19983.gafpack -o na19983_genotype_masked -i na19983 -m cyp2.mask.txt
```

It seems we still can't get the exact genotype, even if the predicted one is close: one haplotype matches the expectations, while the other has low divercency from the true haplotype (how would you calculate this?).

## 6. Graph decomposition

We will use `vg deconstruct` to decompose the graph into individual variants and convert the `gfa` back to a `vcf` file. Additionally to make things easier we will only keep structural variants of at least 100 bp.

```bash
sed 's/#/_/g' cyp2.gfa > cyp2_renamed.gfa
vg deconstruct -P grch38 cyp2_renamed.gfa > cyp2_renamed.vcf
sh filter_variants.sh cyp2_renamed.vcf cyp2_renamed_filt.vcf
```

<details>
<summary> 
<b>filter_variants.sh</b> 
</summary>

```bash 
#!/bin/bash

# Script to filter VCF for variants with at least one allele >= 100bp using bcftools

input_vcf="$1"
output_vcf="$2"

# Filter using bcftools
{
    bcftools view -h "$input_vcf"
    bcftools view -H "$input_vcf" | awk '
    {
        # Check REF allele length
        if (length($4) >= 100) {
            print
            next
        }
        
        # Check ALT allele(s) length
        split($5, alt_alleles, ",")
        for (i in alt_alleles) {
            if (length(alt_alleles[i]) >= 100) {
                print
                next
            }
        }
    }'
} > "$output_vcf"
```
</details>   

Try to have a look at the different alleles. How many alleles does each variant have?

`vg deconstruct` has a parameter called `-L` which allows you to cluster together similar alleles.

```bash
vg deconstruct -P grch38 cyp2_renamed.gfa -L 0.98 > cyp2_renamed_clustered.vcf
sh filter_variants.sh cyp2_renamed_clustered.vcf cyp2_renamed_clustered_filt.vcf
```

How many alleles now?

## 7. Graph construction with vg

We can use our selected variants to build a new graph with `vg` starting from the reference.

```bash
odgi paths -i cyp2_renamed.gfa -L | grep grch38 > ref_path_renamed.txt
odgi paths -K ref_path_renamed.txt -f -o /dev/null -i cyp2_renamed.gfa > ref_path.fa
vg construct -r ref_path.fa -v cyp2_renamed_clustered_filt.vcf -f -A -S -m 100000000 > new.vg
vg view -g new.vg > new.gfa
```

We can use [Bandage](https://rrwick.github.io/Bandage/) to compare the old graph to this new one. What are the differences? (for the first one, put node width to at least 50 and zoom to 1%)
