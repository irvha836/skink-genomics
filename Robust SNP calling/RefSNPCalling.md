# Trimming reads 
```
module load cutadapt
```
```
cutadapt -j 8 -a CCGAGATCGGAAGAGC   -q 25 -o  SQ2252_HMGKKDRX3_s_1_trimmed.fastq.gz --minimum-length 30  SQ2252_HMGKKDRX3_s_1_fastq.txt.gz
cutadapt -j 8 -a CCGAGATCGGAAGAGC   -q 25 -o  SQ2252_HMGKKDRX3_s_2_trimmed.fastq.gz --minimum-length 30  SQ2252_HMGKKDRX3_s_2_fastq.txt.gz 
cutadapt -j 8 -a CCGAGATCGGAAGAGC   -q 25 -o  SQ2252_HMGN2DRX3_s_1_trimmed.fastq.gz --minimum-length 30  SQ2252_HMGN2DRX3_s_1_fastq.txt.gz
cutadapt -j 8 -a CCGAGATCGGAAGAGC   -q 25 -o  SQ2252_HMGN2DRX3_s_2_trimmed.fastq.gz --minimum-length 30  SQ2252_HMGN2DRX3_s_2_fastq.txt.gz
```
```
cd raw
ln -s ../*trimmed* .
```
# Demultiplexing 
```
module load Stacks
process_radtags  -p raw/ -o ./samples/ -b barcodes.txt -e pstI  -c -q --inline_null -r
```

# indexing
```
module load SAMtools
module load BWA 
module load Stacks
```

```
mkdir -p bam

for fq in samples/*.fq.gz
do
    sample=$(basename "$fq" .fq.gz)   # strip extension
    echo "Aligning $sample ..."
    
    bwa mem -t 8 whitakergenome.fasta "$fq" \
    | samtools view -bS - \
    | samtools sort -o bam/${sample}.bam
    
    samtools index bam/${sample}.bam
done
```
# SNP calling 
```
module load Stacks/2.67-GCC-12.3.0
gstacks -I bam/ -M popmap.txt -O gstacks_out/ -t 16
populations -P gstacks_out/ -M popmap.txt -O populations_out/ --vcf
```
Stacks Manual
Julian Catchen1, Nicolas Rochette1, Angel G. Rivera-Colón1,2, William A. Cresko2, Paul A. Hohenlohe3, Angel Amores4, Susan Bassham2, John Postlethwait4

1Department of Evolution, Ecology, and Behavior
University of Illinois at Urbana-Champaign
Urbana, Illinois, 61820
USA
# filtering
```
module load VCFtools
vcftools --vcf populations_out/populations.snps.vcf --missing-indv
```
```
sort -k 5n out.imiss > out.imiss.sorted
```
can rerun this without write single snp and max observe het and R
```
populations -P gstacks_out/ -M popmap.txt \
-O populations_filtered
```
Cut of 90% individuals 


Filtered by pop map, not vcf, Cut of individuals missing 90% of the data but keeping admixed individuals and some founders
(here are individuals i removed)
D01,F09,D09,C07,C06,C05,D29,F02,F28,F13,F10,D04,F25,DN75,PS44

```
mkdir -p populations_filtered90jan22
populations -P gstacks_out/ \
-M popmap90.txt \
-O populations_filtered90jan22/ \
--vcf \
--max-obs-het 0.6 \
-R 0.6
```
Gave me 65451 SNPs!

max missing 0.9 keep
max missing 0.99 (Sup) 











```
mkdir -p populations_filtered2
populations -P gstacks_out/ \
-M popmap.txt \
-O populations_filtered2/ \
--vcf \
--max-obs-het 0.6 \
-R 0.6
```
got me 4403 SNPs
```
populations -P gstacks_out/ -M popmap.txt \
  --write-single-snp --max-observed-het 0.65 -R 0.7 \
  -O populations_filteredv2
```
got me 2064 SNPs

so ludo said try again with R 0.65 

```
mkdir -p populations_filtered2
populations -P gstacks_out/ \
-M popmap.txt \
-O populations_filtered2/ \
--vcf \
--max-obs-het 0.6 \
-R 0.6
```

```
vcftools --vcf populations.snps.vcf \
  --minDP 3 \
  --max-missing 0.6 \
  --remove-indv A04 \
  --remove-indv F10 \
  --remove-indv D04 \
  --remove-indv F25 \
  --remove-indv DN75 \
  --remove-indv PS44 \
  --recode --out robust_filtered_pedigree
```

```
vcftools --vcf populations.snps.vcf \
  --minDP 3 \
  --max-missing 0.6 \
  --remove-indv A04 \
  --remove-indv F10 \
  --remove-indv D04 \
  --remove-indv F25 \
  --remove-indv DN75 \
  --remove-indv PS44 \
  --remove-indv F13 \             
  --recode --out robust_filtered_pedigreev2
```
# phylogeny 
```
git clone https://github.com/edgardomortiz/vcf2phylip.git

# Go into the directory
cd vcf2phylip
```
```
module load python

python3 vcf2phylip.py -i robust_filtered90_minDP3_maxmiss60.recode.vcf -o robust_90filtered_pedigree.phy

```
run IQtree note i also ran without asc model to compare 
```
#!/bin/bash -e
#SBATCH --job-name=iqtree_snp
#SBATCH --output=iqtree_snp.out
#SBATCH --error=iqtree_snp.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Load IQ-TREE module
module load IQ-TREE

# Run IQ-TREE SNP-aware phylogeny
iqtree2 -nt 16 \
  -s robust_filtered90_minDP3_maxmiss60.recode.min4.phy \
  -m GTR+ASC+G \
  -st DNA \
  -bb 1000 \
  -pre robust_filtered90_snp
```


General citation for IQ-TREE 2:

B.Q. Minh, H.A. Schmidt, O. Chernomor, D. Schrempf, M.D. Woodhams, A. von Haeseler, R. Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., 37:1530-1534. https://doi.org/10.1093/molbev/msaa015


# Fresh start for filtering 

Filtered by pop map, not vcf, cut of individuals missing 90% of the data, but keeping admixed individuals and some founders
(here are individuals i removed)
D01,F09,D09,C07,C06,C05,D29,F02,F28,F13,F10,D04,F25,DN75,PS44

```
mkdir -p populations_filtered90jan22
populations -P gstacks_out/ \
-M popmap90.txt \
-O populations_filtered90jan22/ \
--vcf \
--max-obs-het 0.6 \
-R 0.6
```
Gave me 65451 SNPs!

further investigated missingness and it looks good majority of data is below 50% missingness

```
module load VCFtools
vcftools --vcf populations_out/populations.snps.vcf --missing-indv
```
```
sort -k 5n out.imiss > out.imiss.sorted
```
Following SNP calling, genotypes with read depth <3 were removed, and loci were retained only if present in at least 60% of individuals. This ensured consistency between genotype- and locus-level filtering criteria. producing 10460 SNPs (Code below)
```
vcftools --vcf populations.snps.vcf \
  --minDP 3 \
  --max-missing 0.6 \
  --recode \
  --out robust_filtered90_minDP3_maxmiss60

```

I then investigated the relationship between heterozygosity and depth in R within my data set and found no correlation between depth and heterozygosity so we have detected true differences

```
vcftools --vcf robust_filtered90_minDP3_maxmiss60.recode.vcf --het
```
```
vcftools --vcf robust_filtered90_minDP3_maxmiss60.recode.vcf --depth
```
Then re-ran IQ tree (ASK LUDO ABOUT invariant sites) 

```
#!/bin/bash -e
#SBATCH --account=uoo04250
#SBATCH --job-name=iqtree_snp
#SBATCH --output=iqtree_snp.out
#SBATCH --error=iqtree_snp.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Load IQ-TREE module
module load IQ-TREE

# Run IQ-TREE SNP-aware phylogeny
iqtree2 -nt 16 \
  -s robust_filtered90_minDP3_maxmiss60.recode.min4.phy \
  -m GTR+G \
  -st DNA \
  -bb 1000 \
  -pre robust_filtered_snp
```



Annabell stuff 
```
mkdir -p populations_sexesfiltered
populations -P gstacks_out/ \
-M popmap_sex_well_sequenced10F9M.txt \
-O populations_sexesfiltered/ \
--vcf \
--max-obs-het 0.6 \
-R 0.6
```





# FST 
more stringent filtering for FST calculations
```
populations -P gstacks_out/ \
  -M popmap90.txt \
  -O populations_singleSNP/ \
  --vcf \
  --write-single-snp \
  -R 0.6 \
  --max-obs-het 0.6
```
32881 variant sites remained

keep filtering apart from single SNP consistent 
```
vcftools --vcf populations.snps.vcf \
  --minDP 3 \
  --max-missing 0.6 \
  --recode \
  --out robust_filtered90_minDP3_maxmiss60_single
```
kept 4791 SNPs
