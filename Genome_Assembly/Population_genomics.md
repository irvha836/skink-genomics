# PSMC 

need to index my whitakergenome.fasta
and align illumina reads to it with samtools and bwa 

module load SAMtools 
module load BWA 

bwa index whitakergenome.fasta





QC
samtools flagstat whitaker.bam
samtools idxstats whitaker.bam




generate a vcf file from my genome and indexed bam file 

```
module load GCC/11.3.0  # make sure this matches a working version
module load BCFtools/1.16
module load SAMtools/1.16

# Generate VCF from BAM
bcftools mpileup -f whitakergenome.fasta skink_aln.bam | \
bcftools call -c -Ov -o skink.raw.vcf
```
gunzip the large vcf file 

```
module load HTSlib/1.16
bgzip -@ 8 skink.raw.vcf     # creates skink.raw.vcf.gz
tabix -p vcf skink.raw.vcf.gz
```

generate a snp file from the vcf file which will be used downstream for PMSC.
```
module load BCFtools

bcftools view -v snps skink.raw.vcf.gz -Oz -o skink.snps.vcf.gz
tabix -p vcf skink.snps.vcf.gz
```
now filter snps by coverage depth of 57.6 and a minimum of 5

```
bcftools filter -i 'QUAL>20 && DP>20 && DP<300' \
  skink.snps.vcf.gz -Oz -o skink.filtered.vcf.gz
tabix -p vcf skink.filtered.vcf.gz
```

filter by excluding everything that falls outside of the heterozygous ratio e.g outside of a 70/30% site frequences for a or b allele

```
bcftools filter -i 'DP4[2]+DP4[3]>0 && (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>=0.3 && (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])<=0.7' \
  skink.filtered.vcf.gz -Oz -o skink.het.vcf.gz

tabix -p vcf skink.het.vcf.gz
```

/nesi/nobackup/uoo04250/genome_assembly/PSMC/psmc/utils/fq2psmcfa -q20 skink_for_psmc.fq > skin
k.psmcfa
15:47:21 jupyterlab-uoo04250-4cayyfkr /nesi/nobackup/uoo04250/genome_assembly/PSMC $ bcftools consensus -f whitakergenome.fasta skink.het.vcf.gz -H 1 > skink.het.fa
Applied 2367 variants
15:50:57 jupyterlab-uoo04250-4cayyfkr /nesi/nobackup/uoo04250/genome_assembly/PSMC $ /nesi/nobackup/uoo04250/genome_assembly/PSMC/psmc/utils/fq2psmcfa -q20 skink.het.fa > skink.psmcfa
15:56:51 jupyterlab-uoo04250-4cayyfkr /nesi/nobackup/uoo04250/genome_assembly/PSMC $ ls -lh skink.psmcfa
head -n 20 skink.psmcfa
-rw-r--r--. 1 irvha836 uoo04250 14M Sep  2 15:56 skink.psmcfa
>contig_4503_pilon_pilon
NNNNNNNNNTTNNNNNTTTTTTTTTTTNNNTNNNNNNNNNNNTTTTTTTTTTNTNTNTNN
NNNNTTTTTNTNNTNTTNNTNTTNNTTTTTTNNNNNNNNNNNNNNNNNNNNNNNNTNTTN
NNNNNNNTTNNNNNNNNNNNNNNNTTTTTTTTTTTNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNTNNNNNTTTTTTTTTTNNNNNTTTTTTTTNTTTNNNNNTTN
TNNNNTNNNNNNNNNNNNNNNNNNNNNNNTNNNNNNNNTTTTTTTTNTTTTTTNNNTTTT
TNTTTTNTNNNNNNNNNNNTNNTNNTTNTTTTNTNNNNTTTTTTTTTNTTTTTTNTTTTN
NNNNNNNTTTTTTTTTTTTNNNNNNNNTTTNNNNNNNNTTTTTTTTNTTNNNNNTNNTTN
NTTTNTTTNNNNNTTNNNNNNNNNNNNNNNNNNNNNTNTTTTTTTTTTTTTNNNNNNNNT
TTTTTTTTTTNTTTNTNTTNNTTTTTTNNNNTTTTTTTTTNTTTTTTNTTT
>contig_15137_pilon_pilon
TTTTTTTTTTTTTTTTTTTTTTTTTTTTNNNNNNNNNNNTTTTTTTNNNTTTTTNNNNNN
NTTTNTTTTTTTTTTTTNNNTTNTTTTTTTTTNTTTNNTTTTTTTNNNNNNNNNNNTNNN
NTTTTNNNNNNNTTTTTTTTNNNNNNNNNNNNNNNNNTTTTTTTTTNNNNNTTTTNNNNT
TTTTTTTTTTTTTTNNTTTTTTTTTTTTTNNNNNNNTTTTTTTTTTTTTTTTTTNNNTTT
NTNNNTNNNNNTTTTTTNNNNNNNNNNNNTTTTTNNNNNNNNNNNNNNNNNNNT
>contig_17256_pilon_pilon
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTNTTTTTTTNNNNNNNTT
TTTTTNNTNTTTTTTTTTNNNNNNNNTTTTTTNNNNTTTTTTTTTNNNTTTTTTTNNNTT
TTTTNNNNTTTTTTTTTTTTTTTNTTNNTTTTTTTTTTTNNTTTTTTTTTTTTTTTTTTT
16:01:41 jupyterlab-uoo04250-4cayyfkr /nesi/nobackup/uoo04250/genome_assembly/PSMC $ /nesi/nobackup/uoo04250/genome_assembly/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o skink.psmc skink.psmcfa






Next step: build the consensus diploid FASTQ
```
# Make diploid consensus fastq (with IUPAC codes for hets)
bcftools view -v snps skink.filtered.vcf.gz | \
  vcfutils.pl vcf2fq -d 20 -D 300 > skink.fq
```







# psmc trial 2

ran into file conversion errors and lots of missing sites so i used my less masked genome, and reads with higher accuracy (Q20) gave me a coverage of about 27 and a average depth of 22.5


1. mapped reads to genome with bcftools minimap 

```
#!/bin/bash -e
#SBATCH --cpus-per-task  8
#SBATCH --job-name       MapSkink
#SBATCH --mem            64G
#SBATCH --time           72:00:00
#SBATCH --account        uoo04250
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module load SAMtools/1.16.1-GCC-11.3.0
module load minimap2/2.24-GCC-11.3.0

GENOME=pilongapsr2lepido.fasta.masked
READS=reads_q20_l1kv2.fastq.gz

# Align and sort
minimap2 -ax sr -t ${SLURM_CPUS_PER_TASK} $GENOME $READS \
| samtools sort -@ 8 -O BAM -o skink_reads.bam
```

calculate average depth to set up d and D paramaters (it was 21.02) 
```
samtools depth -a skink_reads.bam | \
awk '{sum+=$3} END { print "Average depth = ",sum/NR }'
```

based off this Here option -d sets and minimum read depth and -D sets the maximum. It is
recommended to set -d to a third of the average depth and -D to twice. d should be 7 and D should therefore be 42 and obviously used q scores of 20

so using those paramaters I first indexed the genome, and used the output .bam file to call the variants with bcftools mpileup
Script for bcftools mpileup

```
#!/bin/bash -e
#SBATCH --cpus-per-task  12
#SBATCH --job-name       SkinkPSMC
#SBATCH --mem            32G
#SBATCH --time           24:00:00
#SBATCH --account        uoo04250
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load psmc
module load SAMtools/1.16
module load BCFtools/1.16

GENOME=pilongapsr2lepido.fasta.masked
BAM=skink_reads.bam

# index genome
samtools faidx $GENOME

# index bam
samtools index $BAM

# mpileup -> call -> vcf2fq
bcftools mpileup -Q 20 -q 20 -O v -f $GENOME $BAM \
| bcftools call -c \
| vcfutils.pl vcf2fq -d 7 -D 42 -Q 20 > skink_forpsmc.fq
```

now transform fq file into psmc and perform 100 bootstraps
generation age based on reptiles of nz by hitchmough and van winkle (age of sexual maturity) 6.5 
and per generation 1.17 https://www.nature.com/articles/s41586-023-05752-y#MOESM5 

```
#!/bin/bash -e
#SBATCH --cpus-per-task  6
#SBATCH --job-name       Bt_PSMCHec
#SBATCH --mem            16G
#SBATCH --time           06:00:00
#SBATCH --account        uoo04250
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load psmc

fq2psmcfa skink_forpsmc.fq > skink.psmcfa

splitfa skink.psmcfa > split.psmcfa

psmc -N25 -t15 -r5 -d -p "4+25*2+4+6" -o skink.psmc skink.psmcfa

seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc split2m.psmcfa | sh

cat skink2m.psmc round-*.psmc > skink_bootstrap.psmc

psmc_plot.pl -R -u 1.17e-08 -g 6.5 -p skink_bootstrap skink_bootstrap.psmc
```






