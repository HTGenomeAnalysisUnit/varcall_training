
**Step 2:**  

Retrieve only biallelic variants
```
bcftools view -m2 -M2 -v snps,indels ../../solutions/combined_vcfs/cohort.joint.vcf.gz -Oz -o CASE_2026/cohort.biallelic.vcf.gz

tabix -p vcf cohort.biallelic.vcf.gz
```

Then filter for varaints present in the son but not in either of the parents

```
bcftools view CASE_2026/cohort.biallelic.vcf.gz \
  -i 'GT[0]!="0/0" && GT[1]="0/0" && GT[2]="0/0"' \
  -Oz -o CASE_2026/impossible.raw.vcf.gz

tabix -p vcf CASE_2026/impossible.raw.vcf.gz
```

**Step 3:** 
```
bcftools view -H CASE_2026/impossible.raw.vcf.gz | wc -l
```
**Step 4:**  
```
bcftools view CASE_2026/impossible.raw.vcf.gz \
  -i 'FORMAT/DP[0]>=15 && FORMAT/DP[1]>=15 && FORMAT/DP[2]>=15 && FORMAT/GQ[0]>=30 && FORMAT/GQ[1]>=30 && FORMAT/GQ[2]>=30' \
  -Oz -o CASE_2026/impossible.dp15.gq30.vcf.gz

tabix -p vcf impossible.dp15.gq30.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%AD:%DP:%GQ]\n' \
CASE_2026/impossible.dp15.gq30.vcf.gz > impossible.casefile.tsv
```