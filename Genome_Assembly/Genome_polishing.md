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

Convert fastq Illumina reads to bam files to be aligned to genome to polish.






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
bwa mem -t 16 assembly.fasta forward/merged_forward.fastq.gz reverse/merged_reverse.fastq.gz | samtools view -Sb - > ialigned.bam
```
```
pilon -- assembly.fasta [--frags paired reads.bam] maybe need to do reverse and forward. 

--outdir directory pilon_polished
