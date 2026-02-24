# Admixture 
create admixture environment

```
module load Miniconda3
conda create -n admixture_env python=3.8
#proceed with y
```
activate the environment 
```
conda activate admixture
```
add correct channels
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
install admixture 
```
conda install admixture
```
check it worked (should see admixture there with info about it
```
admixture --help
``` 




generate plink files after removing wild samples as we are only testing for admixed captive individuals
```
module load PLINK/1.09b6.16

plink --vcf robust_filtered90_minDP3_maxmiss60.recode.vcf --make-bed --out robust_filteredadmix --allow-extra-chr
```
prune for linkage disequilibrium
```
plink --
bfile robust_admix_captive \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.2 \
  --out robust_captive_LD
```

```
plink --bfile robust_admix_captive \
  --allow-extra-chr \
  --extract robust_captive_LD.prune.in \
  --make-bed \
  --out robust_captive_pruned
```
retained 3028 SNPs out of 10460.

create file with the proper chromosome format
```
awk '{$1=1; print}' OFS='\t' robust_captive_pruned.bim > fixed.bim
mv fixed.bim robust_captive_pruned.bim
```

Run admixture 
```
for K in $(seq 1 6); do
  admixture --cv=10 -j8 robust_captive_pruned.bed $K | tee robust_captive_K${K}.log
done
```
Check which k model best fits 
```
grep -h "CV error" robust_captive_K*.log
```

CV error (K=1): 0.41866
CV error (K=2): 0.34493
CV error (K=3): 0.31018
CV error (K=4): 0.31119
CV error (K=5): 0.32294
CV error (K=6): 0.32794
therefore k 3 is the best so i should use this however becasue we are working with captive pop from 2 populations  k 2 would be fine


plot in r using robust_captive_pruned.fam and k2 and k3 output.Q files
