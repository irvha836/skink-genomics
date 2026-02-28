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
