Generate a new popmap with harsh filtering (r80 and max-obs-het-65%)for better estimates, aiming for 300-700 SNPs my castle popmap excluded low coverage individuals C05, C06, C07 and included the legacy sample VUW05
```
mkdir -p castle_population_filtered
populations -P gstacks_out/ \
-M castle_popmap.txt \
-O castle_population_filtered/ \
--vcf \
--max-obs-het 0.65 \
-R 0.8
```
Further filter for minDP of 5 and max missing of 0.8 
```
vcftools --vcf populations.snps.vcf \
  --minDP 5 \
  --max-missing 0.8 \
  --recode \
  --out castle_filtered_minDP5_maxmiss80
```
gave me 169 snps 
