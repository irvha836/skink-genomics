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
```
module load Stacks/2.67-GCC-12.3.0
gstacks -I bam/ -M popmap.txt -O gstacks_out/ -t 16
populations -P gstacks_out/ -M popmap.txt -O populations_out/ --vcf
```

```
module load VCFtools
vcftools --vcf populations_out/populations.snps.vcf --missing-indv
```
```
sort -k 4n out.imiss > out.imiss.sorted
```
```
populations -P gstacks_out/ -M popmap.txt \
  --max-obs-het 0.65 -R 0.75 \
  --write-single-snp -O populations_filtered
```
```
vcftools --vcf populations.snps.vcf \
  --minDP 3 \
  --maf 0.005 \
  --max-missing 0.6 \
  --remove-indv A04 \
  --remove-indv F10 \
  --remove-indv D04 \
  --remove-indv F25 \
  --remove-indv DN75 \
  --remove-indv PS44 \
  --recode --out robust_filtered_pedigree
```
