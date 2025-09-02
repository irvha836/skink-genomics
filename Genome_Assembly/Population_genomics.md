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

Next step: build the consensus diploid FASTQ
```
# Make diploid consensus fastq (with IUPAC codes for hets)
bcftools view -v snps skink.filtered.vcf.gz | \
  vcfutils.pl vcf2fq -d 20 -D 300 > skink.fq
```
