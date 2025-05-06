# **SNP Calling**

get to Correct Directory

```
cd uoo04250/popgen 
```

## **Quality Checking** 

check all 4 of the fastq files by generating test small files. (zcat decompresses gz files)
```
zcat SQ2252_HMGKKDRX3_s_1_fastq.txt.gz | head -n 1000000 > test_small_KK_S1.fastq
zcat SQ2252_HMGKKDRX3_s_2_fastq.txt.gz | head -n 1000000 > test_small_KK_S2.fastq
zcat SQ2252_HMGN2DRX3_s_1_fastq.txt.gz | head -n 1000000 > test_small_N2_S1.fastq
zcat SQ2252_HMGN2DRX3_s_2_fastq.txt.gz | head -n 1000000 > test_small_N2_S2.fastq
```

Perform quality check with FastQC 
```
module load FastQC
fastqc test_small_*.fastq
```
## **Trimming adapters and short reads**

cutadapt works by using cutadapt -a ADAPTER [Sequence] [-o output.fastq file] input.fastq file 
use test_small file to look at overrepresented sequences as likely its the adapter will be overrepresented - if you havent got all the adapter then cutadapt will tell you that in majority of samples there is a A T C or G following likely indicating it is part of the adapter
```
module load cutadapt
# Origonally ran this cut adapt code however it didnt work as ive got the wrong adapter in there so obviously it wont do anything.
# Process file 1 (KK S1)
cutadapt -j 8 -a CCGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTCATAGGATCTCG -q 25 --minimum-length 42 --length 42 -o SQ2252_HMGKKDRX3_s_1_trimmed.fastq.gz SQ2252_HMGKKDRX3_s_1_fastq.txt.gz

# Process file 2 (KK S2)
cutadapt -j 8 -a CCGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTCATAGGATCTCG -q 25 --minimum-length 42 --length 42 -o SQ2252_HMGKKDRX3_s_2_trimmed.fastq.gz SQ2252_HMGKKDRX3_s_2_fastq.txt.gz

# Process file 3 (N2 S1)
cutadapt -j 8 -a CCGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTCATAGGATCTCG -q 25 --minimum-length 42 --length 42 -o SQ2252_HMGN2DRX3_s_1_trimmed.fastq.gz SQ2252_HMGN2DRX3_s_1_fastq.txt.gz

# Process file 4 (N2 S2)
cutadapt -j 8 -a CCGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTCATAGGATCTCG -q 25 --minimum-length 42 --length 42 -o SQ2252_HMGN2DRX3_s_2_trimmed.fastq.gz SQ2252_HMGN2DRX3_s_2_fastq.txt.gz
```

# **Demultiplexing**

```
mkdir raw samples
```

```
cd raw 
ln -s ../SQ2252_HMGN2DRX3_s_1_trimmed.fastq.gz .
ln -s ../SQ2252_HMGN2DRX3_s_2_trimmed.fastq.gz .
ln -s ../SQ2252_HMGKKDRX3_s_1_trimmed.fastq.gz .
ln -s ../SQ2252_HMGKKDRX3_s_2_trimmed.fastq.gz .
cd ..
```

 Run `process_radtags` to filter and demultiplex the RADSeq data

```
module load Stacks
process_radtags -p raw/ -o ./samples/ -b barcodes.txt -e pstI -c -q --inline_null -r
```
create popmap file with my barcodes

```
awk '{print $1"\tpop"}' barcodes.txt > popmap.txt
```

had to reformate code above as it used barcodes instead of sample id this fixed it

```
awk 'NR==FNR {map[$1] = $2; next} {print map[$1] "\t" $2}' barcodes.txt popmap.txt > fixed_popmap.txt
```

replace old pop map with fixed one

```
mv fixed_popmap.txt popmap.txt
```
## **Submit snp calling slurm job**

```
for i in 2
do 
    mkdir -p M$i
    echo '#!/bin/sh' > runM$i.sh
    echo "module load Stacks/2.61-gimkl-2022a" >> runM$i.sh
    echo "denovo_map.pl --samples samples/ --popmap popmap.txt -o M$i -M $i -n $i -m 3 -T 16" >> runM$i.sh
    sbatch -A uoo04250 -t 20:00:00 -J M$i -c 16 --mem=64G runM$i.sh
done
```
check job is running

```
squeue -u irvha836 
```
## **Quality check SNP data**

```
module load Stacks
populations -P M2/ -M popmap.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -R 0.5  --write-single-snp -O M2

module load VCFtools 
vcftools --vcf M2/populations.snps.vcf --missing-indv 
sort -k 4n out.imiss
```
shows how much data is missing for each individual 

Identified individuals DN75, F25, PS44, F10 & D04 with over 97% of their data missing, so remove them
can also play around with max missing 

Second run individuals with over 85% of data missing were removed (DN75, PS44, F25, DO4, F10, F21, F28, F02, D29, A03, C05, F23, D09, C06) 

```
vcftools --vcf M2/populations.snps.vcf \
--max-missing 0.6 --minDP 3 \
--remove-indv DN75 --remove-indv F25 --remove-indv PS44 --remove-indv F10 --remove-indv F13 --remove-indv A04 --remove-indv D04 \
--remove-indv F21 --remove-indv F28 --remove-indv F02 --remove-indv D29 --remove-indv A03 --remove-indv C05 --remove-indv F23 \
--remove-indv D09 --remove-indv C06 --maf 0.0001 --recode --out robust_mindp3_r06
```
After filtering and lowering my minDP to 3 from 5 it, kept 3897 out of a possible 3912 Sites
Additonally changed max missing to 0.6 to 0.5 kept 1552 out of a possible 3912 Sites

3795 out of a possible 3912 Sites

testing different minDP paramaters 
```
vcftools --vcf M2/populations.snps.vcf \
--max-missing 0.6 --minDP 4 \
--remove-indv DN75 --remove-indv F25 --remove-indv PS44 --remove-indv F10 --remove-indv F13 --remove-indv A04 --remove-indv D04 \
--remove-indv F21 --remove-indv F28 --remove-indv F02 --remove-indv D29 --remove-indv A03 --remove-indv C05 --remove-indv F23 \
--remove-indv D09 --remove-indv C06 --maf 0.0001 --recode --out robust_mindp4_r06
```
```
vcftools --vcf M2/populations.snps.vcf \
--max-missing 0.6 --minDP 5 \
--remove-indv DN75 --remove-indv F25 --remove-indv PS44 --remove-indv F10 --remove-indv F13 --remove-indv A04 --remove-indv D04 \
--remove-indv F21 --remove-indv F28 --remove-indv F02 --remove-indv D29 --remove-indv A03 --remove-indv C05 --remove-indv F23 \
--remove-indv D09 --remove-indv C06 --maf 0.0001 --recode --out robust_mindp5_r06
```
```
vcftools --vcf M2/populations.snps.vcf \
--max-missing 0.6 --minDP 6 \
--remove-indv DN75 --remove-indv F25 --remove-indv PS44 --remove-indv F10 --remove-indv F13 --remove-indv A04 --remove-indv D04 \
--remove-indv F21 --remove-indv F28 --remove-indv F02 --remove-indv D29 --remove-indv A03 --remove-indv C05 --remove-indv F23 \
--remove-indv D09 --remove-indv C06 --maf 0.0001 --recode --out robust_mindp6_r06
```

### Investigating Depth & Heterozygosity

run independent analysis of Het & Depth on my filtered output file e.g step above run with minDP of 3 and max missing of 0.6 (robust_mindp3_r06.recode.vcf)

```
vcftools --vcf robust_mindp3_r06.recode.vcf --het
```
```
vcftools --vcf robust_mindp3_r06.recode.vcf --depth
```

playing around with different paramaters
```
vcftools --vcf robust_mindp4_r06.recode.vcf --het --out robust_mindp4_r06_het

vcftools --vcf robust_mindp4_r06.recode.vcf --depth --out robust_mindp4_r06_depth
```

```
vcftools --vcf robust_mindp5_r06.recode.vcf --het --out robust_mindp5_r06_het

vcftools --vcf robust_mindp5_r06.recode.vcf --depth --out robust_mindp5_r06_depth
```

```
vcftools --vcf robust_mindp6_r06.recode.vcf --het --out robust_mindp6_r06_het

vcftools --vcf robust_mindp6_r06.recode.vcf --depth --out robust_mindp6_r06_depth
```
