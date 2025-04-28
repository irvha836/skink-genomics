# **SNP Calling**

get to Correct Directory

```
cd popgen 
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

```
module load cutadapt
# Origonally ran this cut adapt code however it didnt work as ive got the wrong adapter in there so obviously it wont do anything.
# Process file 1 (KK S1)
cutadapt -j 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTCATAGGATCTCG -q 25 --minimum-length 55 --length 55 -o SQ2252_HMGKKDRX3_s_1_trimmed.fastq.gz SQ2252_HMGKKDRX3_s_1_fastq.txt.gz

# Process file 2 (KK S2)
cutadapt -j 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTCATAGGATCTCG -q 25 --minimum-length 55 --length 55 -o SQ2252_HMGKKDRX3_s_2_trimmed.fastq.gz SQ2252_HMGKKDRX3_s_2_fastq.txt.gz

# Process file 3 (N2 S1)
cutadapt -j 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTCATAGGATCTCG -q 25 --minimum-length 55 --length 55 -o SQ2252_HMGN2DRX3_s_1_trimmed.fastq.gz SQ2252_HMGN2DRX3_s_1_fastq.txt.gz

# Process file 4 (N2 S2)
cutadapt -j 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTCATAGGATCTCG -q 25 --minimum-length 55 --length 55 -o SQ2252_HMGN2DRX3_s_2_trimmed.fastq.gz SQ2252_HMGN2DRX3_s_2_fastq.txt.gz
```
