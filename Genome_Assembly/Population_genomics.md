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
