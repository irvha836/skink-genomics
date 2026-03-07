Call SNPs from individuals with less than 90% missing data and zero admixture based on admixture, phylogeny and pedigree as can affect estimates 

```
populations -P gstacks_out/ \
-M popmaptrioNorth.txt \
-O TrioML/ \
--vcf \
--max-obs-het 0.65 \
-R 0.8
```
