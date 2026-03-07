Call SNPs from individuals with less than 90% missing data and zero admixture based on admixture, phylogeny and pedigree as can affect estimates 

```
populations -P gstacks_out/ \
-M popmaptrioNorth.txt \
-O TrioML/ \
--vcf \
--max-obs-het 0.65 \
-R 0.7
```
19549 sites remained / 54536 sites 
```
vcftools --vcf populations.snps.vcf \
  --minDP 2 \
  --maf 0.01 \
  --max-missing 0.8 \
  --recode \
  --out northland_filteredtrio_maf001_mm80_mindp2
```
