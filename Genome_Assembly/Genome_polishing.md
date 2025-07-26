**Genome Polishing**
will use Medaka to polish my draft assembly with my ONT reads before polishing further with my illumina reads

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
DRAFT=/nesi/nobackup/uoo04250/genome_assembly/FLYE/FLYEQ8/assembly.fasta
OUTDIR=/nesi/nobackup/uoo04250/genome_assembly/MEDAKA/MedakaONT

# Run Medaka consensus polishing
medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} -m dna_r10.4.1_e8.2_400bps_sup
```


**Polishing with Illumina**

in my illumina data folder 
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
convert paired reads to bam files using bwa aligned to my assembly

```
module load BWA

bwa index assembly.fasta
```
```
module load SAMtools
```

samtools wasnt working so set up my own conda environment 
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
```
conda create -n pilon-env samtools=1.17 bwa=0.7.17 -y
```

```
conda activate pilon-env
```
```
bwa mem -t 16 assembly.fasta trimmed_merged_forward.fastq.gz trimmed_merged_reverse.fastq.gz | samtools view -Sb - > ialigned.bam
```

it was actually 
```
bwa mem -t 16 assembly.fasta trimmed_merged_forward.fastq.gz trimmed_merged_reverse.fastq.gz > aln.sam
```
```
samtools view -bS aln.sam > aln.bam
```
```
samtools sort aln.bam -o aln.sorted.bam
```
```
samtools index aln.sorted.bam
```
pilon -- assembly.fasta [--frags paired reads.bam] maybe need to do reverse and forward. 

--outdir directory pilon_polished


```
#!/bin/bash
#SBATCH --job-name=pilon_polish
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=pilon_%j.out
#SBATCH --error=pilon_%j.err

module load Java/15.0.2
module load Pilon/1.24-Java-15.0.2

java -Xmx16G -jar $EBROOTPILON/pilon.jar \
  --genome assembly.fasta \
  --frags ialigned.sorted.bam \
  --output pilon_round1 \
  --outdir pilon \
  --vcf \
  --changes \
  --threads 16

```


```



#!/bin/bash
#SBATCH --job-name=pilon_round2
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --output=pilon_round2.out
#SBATCH --error=pilon_round2.err

# Load modules
module load Java/1.8.0_144
module load Pilon/1.24

# Create output directory if it doesn't exist
mkdir -p pilon

# Run Pilon (Round 2)
java -Xmx250G -jar $EBROOTPILON/pilon.jar \
  --genome pilon/pilon_round1.fasta \
  --frags ialigned.sorted.bam \
  --output pilon_round2 \
  --outdir pilon \
  --fix snps,indels,gaps \
  --diploid \
  --vcf
```

