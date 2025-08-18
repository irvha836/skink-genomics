# Genome Polishing

## Medaka

need to run medaka which uses my filtered ONT reads to polish my assembly before i polish with illumina reads using pilon

install and activate medika 

```
python3 -m venv medaka_env

pip install --upgrade pip
pip install medaka
source medaka_env/bin/activate
```
```
#!/bin/bash
#SBATCH --job-name=Medaka
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --error=medka_%j.err
#SBATCH --output=medka_%j.out
#SBATCH --time=48:00:00

#activate medaka environment
source /nesi/nobackup/uoo04250/genome_assembly/medaka/medaka_env/bin/activate

# Set the number of threads
NPROC=16

# Set paths to input files and output directory
BASECALLS=/nesi/nobackup/uoo04250/genome_assembly/filtered_reads/all_runs_q8.fastq
DRAFT=/nesi/nobackup/uoo04250/genome_assembly/Illumina/clean_assembly.fasta
OUTDIR=/nesi/nobackup/uoo04250/genome_assembly/MEDAKA/MedakaONT

# Run Medaka consensus polishing
medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} -m dna_r10.4.1_e8.2_400bps_sup
```

# Purge Dups 

need to align my ONT reads to the new consensus fasta alignment produced by medaka
```
#Step 1: Align ONT reads to the consensus genome
minimap2 -x map-ont -t 16 consensus.fasta all_runs_q8_trimmed.fastq > ONTaligned.sam
#Step 2: Convert SAM to BAM
samtools view -Sb ONTaligned.sam > ONTaligned.bam
#Step 3: Sort BAM
samtools sort -@ 16 -o ONTaligned.sorted.bam ONTaligned.bam
#Step 4: Index BAM
samtools index ONTaligned.sorted.bam
```
ran into secondary and suplimentary sequences cloging up purge dups so job wasnt running so had to refilter my sorted alignmed bam to only include primary alignments
```
module load SAMtools
samtools view -F 0x900 -b ONTaligned.sorted.bam > ONTaligned.primary.bam
```
```
module load purge_dups
pbcstat ONTaligned.sorted.bam
```
generate histogram of coverage to purge high coverage duplicate contigs to leave me with only the primary assembly. 
```
module load SAMtools
samtools depth ONTaligned.sorted.bam | awk '{counts[int($3/10)*10]++} END {for (cov in counts) print cov, counts[cov]}' > coverage.hist
```
generate coverage.txt file with 500 10000 and 300000 as my low mid and high coverage depth. 

run purge dups
```
#!/bin/bash
#SBATCH --job-name=purge_dups
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=purge_dups.out
#SBATCH --error=purge_dups.err

# Load required modules
module load minimap2    
module load purge_dups      

# Run pipeline
split_fa consensus.fasta > split.fasta

minimap2 -x asm5 -DP -t 8 split.fasta split.fasta | gzip -c - > self.paf.gz

purge_dups -2 -T cutoffs.txt -c coverage.hist self.paf.gz > dups.bed

get_seqs dups.bed split.fasta > purged.fasta
```
split_fa	Format your genome contigs into one-line FASTA entries
minimap2	Align the genome to itself to detect duplicates
purge_dups	Find and flag likely redundant/haplotig regions
get_seqs	Remove flagged duplicates and output a "purged" genome


**Polishing with Illumina**

prepping illumina reads for polishing (illumina folder) 
```
cat AAH7HYMM5-9684-01-25-01_S1_L001_R1_001.fastq.gz AAH7HYMM5-9684-01-25-02_S2_L001_R1_001.fastq.gz > merged_forward.fastq.gz

cat AAH7HYMM5-9684-01-25-01_S1_L001_R2_001.fastq.gz AAH7HYMM5-9684-01-25-02_S2_L001_R2_001.fastq.gz > merged_reverse.fastq.gz
```
when i ran fastqc it showed reads in tray 1403 were low quality so needed to filter them out/remove them and create a new fastqc 

filtered forward reads 
```
zcat merged_forward.fastq.gz | awk 'BEGIN {n=0} 
/^@.*:1403:/ {n=3; next} 
n > 0 {n--; next} 
{print}' | gzip > filteredf_reads.fastq.gz

#filtered reverse reads 

zcat merged_reverse.fastq.gz | awk 'BEGIN {n=0} 
/^@.*:1403:/ {n=3; next} 
n > 0 {n--; next} 
{print}' | gzip > filtered_reads.fastq.gz
```
then filtered out the random position 6 low quality reads in tray 2202 and 2203 using cut adapt as it was only in the first 6 bp positions so not losing much data  

```
cutadapt -u 0 -U 6 \
  -o trimmed_merged_forward.fastq.gz \
  -p trimmed_merged_reverse.fastq.gz \
  filtered_forward.fastq.gz filtered_reverse.fastq.gz
```

Convert fastq Illumina reads to bam files to be aligned to genome to polish.

identifying weird outlires with bandage

bandage is a software available from https://github.com/rrwick/Bandage/releases

can identify long contigs and isolate sequences for blast from my gfa assembly graph. 
identified one contig as Serratia liquefaciens. so removed using seqkit 

```
module load Seqkit
# Remove by node name
seqkit grep -v -n -p contig_4650 assembly.fasta > clean_assembly.fasta
```
edge_14531 the other long outlier all came back with reptile hits from blast so decided not to remove. 
Anolis and Heteronotia) and a chelonian (Stigmochelys).

checked i removed contig 4650 
```
grep "contig_4640" clean_assembly.fasta
```


**Pilon**
convert paired illumina reads to bam files using bwa aligned to my assembly (gapclosed.fasta)

```
#!/bin/bash -e
#SBATCH --job-name=align_gapclosed
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load modules
module load BWA/0.7.17-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

# Input files
GENOME=gapclosed.fasta
READS1=trimmed_merged_forward.fastq.gz
READS2=trimmed_merged_reverse.fastq.gz
OUTBAM=aln_gapclosed.sorted.bam

# Index genome for BWA
bwa index $GENOME

# Align and convert to sorted BAM
bwa mem -t $SLURM_CPUS_PER_TASK $GENOME $READS1 $READS2 | \
    samtools sort -@ $SLURM_CPUS_PER_TASK -o $OUTBAM -

# Index BAM
samtools index $OUTBAM
```
now run pilon
```
#!/bin/bash
#SBATCH --job-name=pilon_round2
#SBATCH --output=pilon_round2v2.out
#SBATCH --error=pilon_round2v2.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G

# Define paths
PILON_JAR=~/pilon_local/bin/pilon.jar
GENOME=gapclosed.fasta
BAM=aln_gapclosed.sorted.bam
OUT_PREFIX=pilongaps

# Run Pilon
java -Xmx400G -jar $PILON_JAR \
  --genome $GENOME \
  --frags $BAM \
  --output $OUT_PREFIX \
  --vcf \
  --changes \
  --diploid
```
Then run a second time for further polishing

```
#!/bin/bash
#SBATCH --job-name=pilon_round2
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=512G
#SBATCH --output=pilon_round2.out
#SBATCH --error=pilon_round2.err

# Load modules
module load Java/1.8.0_144
module load Pilon/1.24

# Run Pilon (Round 2)
java -Xmx480G -jar $EBROOTPILON/pilon.jar \
  --genome pilon_round1.fasta \
  --frags r1aln.sorted.bam \
  --output pilon_round2 \
  --outdir . \
  --fix snps,indels,gaps \
  --diploid \
  --vcf
```




















had to split polishing for round 2 in half due to memory issues 
```
module load SAMtools
samtools faidx pilon_round1.fasta
cut -f1 pilon_round1.fasta.fai > all_scaffolds.txt
split -n l/2 all_scaffolds.txt scaffolds_
awk '{print $1"\t0\t1000000000"}' scaffolds_aa > targets1.bed
awk '{print $1"\t0\t1000000000"}' scaffolds_ab > targets2.bed
```









#!/bin/bash
#SBATCH --job-name=pilon_round2
#SBATCH --output=pilon_round2v2b.out
#SBATCH --error=pilon_round2v2b.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G

PILON_JAR=~/pilon_local/bin/pilon.jar
OUT_PREFIX=pilon_round2b

cd /nesi/nobackup/uoo04250/genome_assembly/Illumina/pilon

java -Xmx210G -jar $PILON_JAR \
  --genome pilon_round1.fasta \
  --frags r1aln.sorted.bam \
  --output $OUT_PREFIX \
  --vcf \
  --changes \
  --diploid \
  --targets targets2.bed




**BUSCO POLISH**
```
#!/bin/bash
#SBATCH --job-name=busco_compare
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=busco_compare_%j.out
#SBATCH --error=busco_compare_%j.err

cd /home/irvha836/uoo04250/genome_assembly/Illumina/pilon

# Load BUSCO via Apptainer
apptainer exec ./busco_5.8.2--pyhdfd78af_0.sif \
busco -i pilon_round1.fasta \
-l sauropsida_odb10 \
-o busco_round1 \
-m genome \
--cpu 16 -f
```
```
#!/bin/bash
#SBATCH --job-name=busco_compare
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=busco_compare_%j.out
#SBATCH --error=busco_compare_%j.err

cd /home/irvha836/uoo04250/genome_assembly/Illumina/pilon
apptainer exec ./busco_5.8.2--pyhdfd78af_0.sif \
busco -i pilon_round2.fasta \
-l sauropsida_odb10 \
-o busco_round2 \
-m genome \
--cpu 16 -f
```





**Double Check of Completness using merqury** 

#!/bin/bash
#SBATCH --job-name=merqury_run
#SBATCH --output=merqury_run.out
#SBATCH --error=merqury_run.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

module load Merqury/1.3-Miniconda3

# Inputs
trimmed_merged_forward.fastq.gz
trimmed_merged_reverse.fastq.gz
ASSEMBLY=clean_assembly.fasta
K=21

# Outputs
READ_KMER=read_k${K}.meryl
ASM_KMER=asm_k${K}.meryl
MERQURY_OUT=merqury_output

# Step 1: Make k-mers from reads
meryl count k=$K output $READ_KMER $READ_R1 $READ_R2

# Step 2: Make k-mers from assembly
meryl count k=$K output $ASM_KMER $ASSEMBLY

# Step 3: Run Merqury
merqury.sh $READ_KMER $ASSEMBLY $MERQURY_OUT

