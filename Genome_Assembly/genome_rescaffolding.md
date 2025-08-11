# Gap closing with Cobbler and Rails

```
#!/bin/bash -e
#SBATCH --job-name=RailsCobbler
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --account=uoo02423
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --hint=nomultithread

# ====== Module setup ======
module purge
module load Perl/5.34.1-GCC-11.3.0
module load minimap2/2.24-GCC-11.3.0
module load SAMtools/1.13-GCC-9.2.0

# ====== Paths ======
ASSEMBLY="purged_5kb.fa"                  # Assembly FASTA file
LONGREADS="all_runs_q8_trimmed.fastq.gz"  # Long reads file
ANCHOR=250                                # Anchoring length
IDENTITY=0.85                             # Minimum identity
GRACE=500                                 # Grace length
THREADS=12                                # Threads to use
READTYPE="ont"                            # Long read type: ont / map-pb / etc.
SAMTOOLS_PATH="$(which samtools)"         # Auto-detect samtools path

# Full path to local RAILS/Cobbler executables folder
RAILS_DIR="/nesi/nobackup/uoo04250/genome_assembly/gapclosing/RAILS/bin"

# ====== Step 1: Reformat FASTA ======
echo "Reformatting assembly FASTA..."
perl -ne 'if(/^>/){print $_} else {chomp; print "$_\n"}' "$ASSEMBLY" > "${ASSEMBLY%.fa*}-formatted.fa"

# ====== Step 2: Create FOF for long reads ======
echo "$LONGREADS" > longreads.fof

# ====== Step 3: Cobbler ======
echo "Running Cobbler..."
minimap2 -x map-"$READTYPE" -t "$THREADS" "${ASSEMBLY%.fa*}-formatted.fa" "$LONGREADS" \
    | "$SAMTOOLS_PATH" view -bS - \
    | perl "$RAILS_DIR/cobbler.pl" \
        -f "${ASSEMBLY%.fa*}-formatted.fa" \
        -s stream \
        -d "$ANCHOR" \
        -i "$IDENTITY" \
        -g "$GRACE" \
        -q longreads.fof \
        -p "$SAMTOOLS_PATH"

# ====== Step 4: RAILS ======
echo "Running RAILS..."
minimap2 -x map-"$READTYPE" -t "$THREADS" "${ASSEMBLY%.fa*}-formatted.fa" "$LONGREADS" \
    | "$SAMTOOLS_PATH" view -bS - \
    | perl "$RAILS_DIR/RAILS" \
        -f "${ASSEMBLY%.fa*}-formatted.fa" \
        -s stream \
        -d "$ANCHOR" \
        -i "$IDENTITY" \
        -g "$GRACE" \
        -q longreads.fof \
        -p "$SAMTOOLS_PATH"

echo "RAILS+Cobbler minimap2 streaming run complete."

```

# lrgapcloser 

download and extract the tar files from the gapcloser github
```
wget https://github.com/CAFS-bioinformatics/LR_Gapcloser/archive/refs/heads/master.tar.gz -O LR_Gapcloser.tar.gz
tar -zxvf LR_Gapcloser.tar.gz
```

#!/bin/bash -e
#SBATCH --cpus-per-task=12
#SBATCH --job-name=lr_gapc
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --hint=nomultithread

module load BWA/0.7.17-GCC-11.3.0

export PATH=$PATH:/nesi/nobackup/uoo04250/genome_assembly/gapclosing/lrcloser/LR_Gapcloser-master/src

sh LR_Gapcloser.sh \
  -i /path/to/scaffolds_with_gaps.fasta \
  -l /path/to/long_reads_corrected.fasta \
  -s n \
  -t ${SLURM_CPUS_PER_TASK} \
  -r 10
