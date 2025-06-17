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

then combine all fastq pass files from each run into one directory (concatenating the files) 

```
cat /home/irvha836/uoo04250/genome_assembly/dorado_basecalls/run*/pass/*.fastq > genome_all_pass.fastq
```
