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
