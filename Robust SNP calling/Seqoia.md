Generate a new popmap with harsh filtering (r80 and max-obs-het-65%)for better estimates, aiming for 300-700 SNPs my castle popmap excluded low coverage individuals C05, C06, C07 and included the legacy sample VUW05

# Castle individuals
```
mkdir -p castle_population_filtered
populations -P gstacks_out/ \
-M castle_popmap.txt \
-O castle_population_filtered/ \
--vcf \
--max-obs-het 0.65 \
-R 0.8
```
```
module load VCFtools
vcftools --vcf populations_out/populations.snps.vcf --missing-indv
sort -k 5n out.imiss > out.imiss.sorted
```

Further filter for minDP of 5 and max missing of 0.8 
```
vcftools --vcf populations.snps.vcf \
  --minDP 5 \
  --max-missing 0.8 \
  --recode \
  --out castle_filtered_minDP5_maxmiss80
```
gave me 169 snps out of 4110 sites

lowering mindp to 3
```
vcftools --vcf populations.snps.vcf \
  --minDP 3 \
  --max-missing 0.8 \
  --recode \
  --out castle_filtered_minDP3_maxmiss80
```
gave me 1111 out of a possible 4110 Sites

mindp4 kept 474 out of a possible 4110 Sites

to get snp counts 
```
bcftools view -H castle_filtered_minDP3_maxmiss80.recode.vcf | wc -l
```

# Northland individuals

removed castle and admixed individuals from popmap as well as individuals with more than 90% missing data apart from pedigree important individuals such as F21 & F23 
(D01, F09, B01, A02, D09, D29, F02, F28, A03, F13, F10, A04, D04, F25, DN75, PS44)

```
populations -P gstacks_out/ \
-M northland_popmap.txt \
-O northland_population_filtered/ \
--vcf \
--max-obs-het 0.65 \
-R 0.8
```
```
vcftools --vcf populations.snps.vcf \
  --minDP 5 \
  --max-missing 0.8 \
  --recode \
  --out northland_filtered_minDP5_maxmiss80
 ```

Keeping 57 out of 57 individuals. However, this only retained 5 out of 6898 snps.


Tried again with Mindp3 max-missing 0.8 kept 273 out of a possible 6898 Sites, mindp 2 kept 1283 out of a possible 6898 Sites.
