# Contamination checking with Diamond and uniprot 

```
module purge
module load DIAMOND/2.1.14-GCC-12.3.0

# Build DIAMOND database from Swiss-Prot FASTA
diamond makedb \
  --in uniprot_sprot.fasta \
  -d uniprot_sprot
```


```
#!/bin/bash -e
#SBATCH --job-name=whitakers_diamond
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --account=uoo04250
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --hint=nomultithread

module purge
module load DIAMOND/2.1.14-GCC-12.3.0

DB=/nesi/nobackup/uoo04250/ref_snps/contamination/uniprot_sprot.dmnd
QUERY=/nesi/nobackup/uoo04250/ref_snps/whitakergenome.fasta

diamond blastx \
  --query "$QUERY" \
  --db "$DB" \
  --out whitakers_vs_sprot.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --sensitive \
  --max-target-seqs 1 \
  --evalue 1e-25 \
  --query-cover 10 \
  --threads $SLURM_CPUS_PER_TASK \
  --tmpdir $TMPDIR
```


# Q20 Coverage 
reads were first filtered using nanofilt for q20 and aligned to the genome 

```
#!/bin/bash -e
#SBATCH --job-name=whitakers_ont_map
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --account=uoo04250
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --hint=nomultithread

module purge
module load minimap2/2.26-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0

GENOME=/nesi/nobackup/uoo04250/ref_snps/whitakergenome.fasta
READS=/nesi/nobackup/uoo04250/genome_assembly/filtered_reads/reads_q20_l1kv2.fastq.gz

OUTDIR=/nesi/nobackup/uoo04250/ref_snps/coverage_q20
mkdir -p $OUTDIR

BAM=$OUTDIR/whitakers_q20_mapped.bam

# Map ONT reads
minimap2 -ax map-ont -t $SLURM_CPUS_PER_TASK $GENOME $READS | \
samtools sort -@ $SLURM_CPUS_PER_TASK -o $BAM

# Index BAM
samtools index $BAM

# Mapping statistics
samtools flagstat $BAM > $OUTDIR/mapping_stats.txt

# Depth file (per-base depth)
samtools depth -aa $BAM > $OUTDIR/depth_per_base.txt
```

calculate depth in 100kb windows

```
#!/bin/bash -e
#SBATCH --job-name=whit_cov100kb
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --account=uoo04250
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --hint=nomultithread

module purge
module load SAMtools/1.19-GCC-12.3.0
module load BEDTools/2.31.1-GCC-12.3.0

GENOME=/nesi/nobackup/uoo04250/ref_snps/whitakergenome.fasta
BAM=/nesi/nobackup/uoo04250/ref_snps/coverage_q20/whitakers_q20_mapped.bam
OUTDIR=/nesi/nobackup/uoo04250/ref_snps/coverage_q20

W=100000

mkdir -p "$OUTDIR"

# Ensure fasta index exists
if [ ! -s "${GENOME}.fai" ]; then
  samtools faidx "$GENOME"
fi

# Genome sizes file
cut -f1,2 "${GENOME}.fai" > "$OUTDIR/genome.sizes"

# Make windows
bedtools makewindows -g "$OUTDIR/genome.sizes" -w $W > "$OUTDIR/windows_${W}.bed"

# Mean depth per window (last column is mean depth)
bedtools coverage -a "$OUTDIR/windows_${W}.bed" -b "$BAM" -mean > "$OUTDIR/coverage_${W}_mean.txt"

# Summary stats
echo "### Coverage summary (window=${W} bp)" > "$OUTDIR/coverage_${W}_summary.txt"

awk '{sum+=$NF; n++} END {print "mean_depth\t" sum/n "\nwindows\t" n}' \
  "$OUTDIR/coverage_${W}_mean.txt" >> "$OUTDIR/coverage_${W}_summary.txt"

awk '{print $NF}' "$OUTDIR/coverage_${W}_mean.txt" | sort -n | \
  awk '{a[NR]=$1} END {print "median_depth\t" a[int((NR+1)/2)]}' \
  >> "$OUTDIR/coverage_${W}_summary.txt"

awk '$NF < 1 {c++} END {print "windows_mean_depth_lt_1x\t" c}' \
  "$OUTDIR/coverage_${W}_mean.txt" >> "$OUTDIR/coverage_${W}_summary.txt"

awk '$NF == 0 {c++} END {print "windows_zero_depth\t" c}' \
  "$OUTDIR/coverage_${W}_mean.txt" >> "$OUTDIR/coverage_${W}_summary.txt"

echo "Wrote:"
echo " - $OUTDIR/coverage_${W}_mean.txt"
echo " - $OUTDIR/coverage_${W}_summary.txt"
```

```
#!/bin/bash -e
#SBATCH --job-name=whit_cov100kb
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --account=uoo04250
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --hint=nomultithread

module purge
module load mosdepth/0.3.4-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

BAM=/nesi/nobackup/uoo04250/ref_snps/coverage_q20/whitakers_q20_mapped.bam
OUTDIR=/nesi/nobackup/uoo04250/ref_snps/coverage_q20/mosdepth_100kb
PREFIX=$OUTDIR/whitakers_q20

mkdir -p "$OUTDIR"

samtools index -@ $SLURM_CPUS_PER_TASK "$BAM"

mosdepth \
  --threads $SLURM_CPUS_PER_TASK \
  --by 100000 \
  --no-per-base \
  --flag 2304 \
  "$PREFIX" "$BAM"
```


# mitogenome coverage 

```
#!/bin/bash -e
#SBATCH --job-name=mito_cov
#SBATCH --account=uoo04250
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=mito_cov_%j.out
#SBATCH --error=mito_cov_%j.err

module purge
module load minimap2
module load SAMtools
module load mosdepth

# ---- input files
REF="whitaker_mitogenomenew.fasta"
READS="reads_q20_l1kv2.fastq.gz"
PREFIX="whitaker_mitonewcoverage"

# ---- map ONT reads to mitogenome
minimap2 -ax map-ont -t ${SLURM_CPUS_PER_TASK} "$REF" "$READS" | \
    samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${PREFIX}.sorted.bam

# ---- index bam
samtools index ${PREFIX}.sorted.bam

# ---- per-base depth with samtools
samtools depth -a ${PREFIX}.sorted.bam > ${PREFIX}.depth.txt

# ---- summary stats from depth file
awk '{sum+=$3; n++} END {print "Mean_coverage\t" sum/n}' ${PREFIX}.depth.txt > ${PREFIX}.coverage_stats.txt

awk 'NR==1{min=$3; max=$3}
     {if($3<min) min=$3; if($3>max) max=$3}
     END {print "Min_coverage\t" min "\nMax_coverage\t" max}' \
     ${PREFIX}.depth.txt >> ${PREFIX}.coverage_stats.txt

awk '$3>=1{a++} $3>=10{b++} $3>=20{c++} {n++}
     END {
       print "Pct_bases_>=1x\t" 100*a/n
       print "Pct_bases_>=10x\t" 100*b/n
       print "Pct_bases_>=20x\t" 100*c/n
     }' ${PREFIX}.depth.txt >> ${PREFIX}.coverage_stats.txt

# ---- mosdepth output too
mosdepth -t ${SLURM_CPUS_PER_TASK} ${PREFIX} ${PREFIX}.sorted.bam
```
