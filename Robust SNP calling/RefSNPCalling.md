# Trimming reads 

module load cutadapt 

cutadapt -j 8 -a CCGAGATCGGAAGAGC   -q 25 -o  SQ2252_HMGKKDRX3_s_1_trimmed.fastq.gz --minimum-length 30  SQ2252_HMGKKDRX3_s_1_fastq.txt.gz
cutadapt -j 8 -a CCGAGATCGGAAGAGC   -q 25 -o  SQ2252_HMGKKDRX3_s_2_trimmed.fastq.gz --minimum-length 30  SQ2252_HMGKKDRX3_s_2_fastq.txt.gz 
cutadapt -j 8 -a CCGAGATCGGAAGAGC   -q 25 -o  SQ2252_HMGN2DRX3_s_1_trimmed.fastq.gz --minimum-length 30  SQ2252_HMGN2DRX3_s_1_fastq.txt.gz
cutadapt -j 8 -a CCGAGATCGGAAGAGC   -q 25 -o  SQ2252_HMGN2DRX3_s_2_trimmed.fastq.gz --minimum-length 30  SQ2252_HMGN2DRX3_s_2_fastq.txt.gz
