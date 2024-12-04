# Structural variant calling

Here we perform structural variant (SV) calling from long-read reads from Pacific Biosciences using [sniffles](https://github.com/fritzsedlazeck/Sniffles)

## Input data

Input folder: `/project/varcall_training/data/partial/alignment/`
Files: `hg19.svs.bam`, `hg19.svs.bam.csi`

## Suggested computational resources

To process the small test dataset, we suggest the following computational resources:

- CPUs: 1
- Memory: 1 GB

## Commands

Navigate to your desired output folder and run the following command:

```bash
singularity exec -B $PWD -B /project/varcall_training -B /localscratch \
	/project/varcall_training/bin/sniffles2_v2.2.sif \
	sniffles \
	-i /project/varcall_training/data/partial/alignment/hg19.svs.bam \
	-v svs.vcf.gz
```

### Options explained

| Option | Description |
|--------|-------------|	
| `-i` | Input .bam |
| `-o` | Output vcf (.gz) and its index |

### Some statistics

Check a couple of lines from the generated vcf file in the output folder with `bcftools view`
