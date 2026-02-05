
```bash
pggb -c 2 -k 101 -x auto -i haplotypes_plus_ref.fa -n 9 -o cyp2_graph
odgi paths -i cyp2_graph/*smooth.final.og -L | grep grch > ref_path.txt  
odgi sort -Y -H ref_path.txt  -i cyp2_graph/*smooth.final.og -o cyp2_graph/cyp.og
odgi view -i cyp2_graph/cyp.og -g > cyp2.gfa
zgrep '^S' cyp2.gfa | awk '{{print("node."$2,length($3))}}' OFS="\\t" > cyp2.node_length.tsv
odgi paths -i cyp2.gfa -H | cut -f 1,4- | gzip > cyp2.node_cov.tsv.gz
bwa-mem2 index haplotypes_plus_ref.fa
bwa-mem2 mem -t 2 -p -h 1000 haplotypes_plus_ref.fa na19983.fa.gz | samtools sort -T tmp | samtools view -o na19983.realigned.bam
gfainject --gfa cyp2.gfa --bam na19983.realigned.bam --alt-hits 10000 | gzip > na19983.gaf
gafpack --gfa cyp2.gfa --gaf na19983.gaf --len-scale --weight-queries | gzip > na19983.gafpack
cosigt -p cyp2.node_cov.tsv.gz -g na19983.gafpack -o na19983_genotype -i na19983
panplexity --input-gfa cyp2.gfa -t auto -k 16 -w 100 -d 100 --complexity linguistic -m cyp2.mask.txt --threads 2
cosigt -p cyp2.node_cov.tsv.gz -g na19983.gafpack -o na19983_genotype_masked -i na19983 -m cyp2.mask.txt
```
