**Genome Aseembly** 

just house keeping making sure all my pod5 files are in easy accessible directories for rebasecalling by copying all pod5 runs into a pod5 folder
e.g run one

```
cp -r /home/irvha836/uoo04250/040325_WSKINK_1/WSKINK_1/20250304_1437_P2S-00650-B_PAY16484_e8780dd9/pod5 \
>       /home/irvha836/uoo04250/genome_assembly/pod5/run1_pod5
```
```
cp -r /home/irvha836/uoo04250/070325_WSKINK_1B/WSKINK_1B/20250307_1554_P2S-00650-B_PAY16484_9d98130a/pod5 \
/home/irvha836/uoo04250/genome_assembly/pod5/run2_pod5

```
make directory for my rebasecalling of my pod5files using dorado

```
mkdir -p /home/irvha836/uoo04250/genome_assembly/dorado_basecalls/
```
```
nano dorado_array_basecall.sh
```

SLURM JOB FOR REBASECALLING. 
```
#!/bin/bash
#SBATCH --job-name=dorado_basecall
#SBATCH --account uoo04250
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --array=1-4
#SBATCH --output=dorado_%A_%a.out

# Load Dorado module
module load Dorado/0.9.1

# Define input, output, and model directories
INPUT_DIR="/home/irvha836/uoo04250/genome_assembly/pod5"
OUTPUT_DIR="/home/irvha836/uoo04250/genome_assembly/dorado_basecalls"
MODEL_DIR="/home/irvha836/uoo04250/genome_assembly/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Map SLURM_ARRAY_TASK_ID to pod5 run subdirectories
case $SLURM_ARRAY_TASK_ID in
    1)
      RUN_DIR="run1_pod5"
      ;;
    2)
      RUN_DIR="run2_pod5"
      ;;
    3)
      RUN_DIR="run3_pod5"
      ;;
    4)
      RUN_DIR="run4_pod5"
      ;;
esac

echo "Starting basecalling for $RUN_DIR"

# Create output directory for this run
mkdir -p $OUTPUT_DIR/$RUN_DIR

# Run Dorado basecalling using local model directory
dorado basecaller "$MODEL_DIR" \
    $INPUT_DIR/$RUN_DIR \
    --recursive \
    --device cuda:all \
    --emit-fastq \
    --output-dir $OUTPUT_DIR/$RUN_DIR

echo "Finished basecalling for $RUN_DIR"
```

Submit SLURM array 

```
sbatch dorado_array_basecall.sh
```
Quality checking run 2 - can either do this individually or independently 

```
module load NanoPlot
```

```
NanoPlot --fastq calls_2025-06-16_T03-25-09.fastq \
         --loglength \
         --N50 \
         --outdir run2_nanoplot_qc \
         --threads 8

```
```
NanoPlot --fastq calls_2025-06-16_T09-29-44.fastq \
         --loglength \
         --N50 \
         --outdir run4_nanoplot_qc \
         --threads 8
```
```
NanoPlot --fastq calls_2025-06-16_T03-25-10.fastq \
         --loglength \
         --N50 \
         --outdir run1_nanoplot_qc \
         --threads 8
```
```
NanoPlot --fastq calls_2025-06-16_T03-25-10.fastq \
         --loglength \
         --N50 \
         --outdir run3_nanoplot_qc \
         --threads 8
```

filter rebasecalled reads by Q score of 8 using Nanofilt
```
module load nanofilt
```

submit as job 
```
#!/bin/bash
#SBATCH --job-name=nanofilt_q8
#SBATCH --output=nanofilt_q8_%A_%a.out
#SBATCH --error=nanofilt_q8_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Load the NanoFilt module (optional, in case dependencies are needed)
module load nanofilt/2.6.0-gimkl-2020a-Python-3.8.2

# Output directory for filtered reads
OUTDIR="/home/irvha836/uoo04250/genome_assembly/filtered_reads"
mkdir -p "$OUTDIR"

# Array of input files
FILES=(/home/irvha836/uoo04250/genome_assembly/dorado_basecalls/run*_pod5/*.fastq)

# Select file based on SLURM array task ID
f=${FILES[$SLURM_ARRAY_TASK_ID]}

echo "Filtering $f for Q>=8"

# Get base filename without extension
base=$(basename "$f" .fastq)

# Set output file path
out="${OUTDIR}/${base}_q8.fastq"

# Use NanoFilt with full path to filter reads with Q >= 8
cat "$f" | /opt/nesi/CS400_centos7_bdw/nanofilt/2.6.0-gimkl-2020a-Python-3.8.2/bin/NanoFilt -q 8 > "$out"

echo "Filtered reads saved to $out"
```
```
sbatch --array=0-3 nanofilt_filter.sh
```


then combine all fastq pass files from each run into one (concatenating the files) 

```
cat *_q8.fastq > all_runs_q8.fastq
```



**FLYE Assembly**

```
#!/bin/bash
#SBATCH --job-name=FLYEQ8
#SBATCH --output=Flyeq8_%A_%a.out
#SBATCH --error=Flyeq8_%A_%a.err
#SBATCH --time=120:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=16
#SBATCH --partition=genoa

module load Flye/2.9.5-foss-2023a-Python-3.11.6

mkdir -p /home/irvha836/uoo04250/genome_assembly/FLYE/FLYEQ8

flye --nano-raw /home/irvha836/uoo04250/genome_assembly/filtered_reads/all_runs_q8.fastq \
     --out-dir /home/irvha836/uoo04250/genome_assembly/FLYE/FLYEQ8 \
     --genome-size 1.6g \
     --threads ${SLURM_CPUS_PER_TASK} \
     --iterations 3
```

```
tail -f flye_assembly_*.out
```
preliminary stats 
```
awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' assembly.fasta | \
sort -nr | awk '{sum+=$1; a[NR]=$1} END {for (i=1;i<=NR;i++) {s+=a[i]; if (s>=sum/2) {print "N50 = " a[i]; exit}}}'
```
N50 = 3635236

```
module load SeqKit
seqkit stats assembly.fasta
```

file            format  type  num_seqs        sum_len  min_len    avg_len     max_len

assembly.fasta  FASTA   DNA      6,602  1,500,389,369      147  227,262.9  25,721,214

**BUSCO**
trying to laod in busco environment
```
apptainer pull docker://quay.io/biocontainers/busco:5.8.2--pyhdfd78af_0
```
```
#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=busco_%j.out
#SBATCH --error=busco_%j.err
cd /nesi/nobackup/uoo04250/genome_assembly/FLYE/FLYEQ8
apptainer exec ./busco_5.8.2--pyhdfd78af_0.sif \
busco -i assembly.fasta \
-l sauropsida_odb10 \
-o busco_output \
-m genome \
--cpu 16 -f
```

BUSCO stats
	C:97.1%[S:93.8%,D:3.2%],F:0.7%,M:2.2%,n:7480,E:2.6%	   
	7260	Complete BUSCOs (C)	(of which 186 contain internal stop codons)		   
	7018	Complete and single-copy BUSCOs (S)	   
	242	Complete and duplicated BUSCOs (D)	   
	52	Fragmented BUSCOs (F)			   
	168	Missing BUSCOs (M)			   
	7480	Total BUSCO groups searched		   

Assembly Statistics:
	6602	Number of scaffolds
	6602	Number of contigs
	1500389369	Total length
	0.000%	Percent gaps
	3 MB	Scaffold N50
	3 MB	Contigs N50
