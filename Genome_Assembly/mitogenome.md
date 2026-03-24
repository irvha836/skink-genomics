1. extract contigs of interest by blasting oligosoma alani mitochondrial genes from NBCI against my genome

```
#!/bin/bash -e
#SBATCH --job-name=mitoBLAST
#SBATCH --time=02:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8
#SBATCH --account=uoo04250

module load BLAST/2.16.0-GCC-12.3.0
module load SAMtools/1.16.1-GCC-11.3.0   # load samtools too

# Step 1: make your genome into a BLAST database
makeblastdb -in whitakergenome.fasta -dbtype nucl -out whitakerDB

# Step 2: run BLASTn of 12S query against genome
blastn \
  -query robust12s.fasta \
  -db whitakerDB \
  -out mito_vs_whitaker.tsv \
  -evalue 1e-10 \
  -num_threads $SLURM_CPUS_PER_TASK \
  -outfmt "6 qseqid qstart qend sseqid sstart send length pident evalue bitscore"

# Step 3: collect contig IDs with hits
cut -f4 mito_vs_whitaker.tsv | sort -u > hit_contigs.txt

# Step 4: extract those contigs from the genome
samtools faidx whitakergenome.fasta $(cat hit_contigs.txt) > candidate_mito_contigs.fasta

echo "Done! Hits saved in mito_vs_whitaker.tsv and sequences in candidate_mito_contigs.fasta"
```

narrowed it down to one contig (32kb) and then looked at where the hits for multiple genes were to see where in the contig the mitogenome was 
```
#!/bin/bash -e
#SBATCH --job-name=mitoGeneBLAST
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=uoo04250

# Load BLAST module
module load BLAST/2.16.0-GCC-12.3.0

# Step 1: make your candidate contigs a BLAST database
makeblastdb -in candidate_mito_contigs.fasta -dbtype nucl -out candidateDB

# Step 2: run BLASTn of your multi-gene mitogenome query
blastn \
  -query alani_mito_genes.fasta \
  -db candidateDB \
  -out mitoGenes_vs_candidate.out \
  -evalue 1e-10 \
  -num_threads $SLURM_CPUS_PER_TASK \
  -outfmt "6 qseqid qstart qend sseqid sstart send length pident evalue bitscore"
```

then looked at coverage and plotted coverage as expected there to be higher coverage due to their being more mitochondria with samtools
```
samtools depth aln.bam | awk '{sum+=$3} END { print "Average depth = ",sum/NR }'
```
```
depth <- read.table("depth.txt")
plot(depth$V2, depth$V3, type="l", xlab="Position", ylab="Depth", main="Coverage across contig_17657")
abline(h=260, col="red", lty=2) # average depth
```
<img width="1436" height="952" alt="image" src="https://github.com/user-attachments/assets/f52edef3-f7f4-470e-8b1f-9ca446766f80" />
then removed start and end of contig 7000bp and 30000bp tail leaving me with 23kb candidate region of high coverage

```
samtools faidx whitakergenome.fasta contig_17657_pilon_pilon:7000-30000 > whitaker_mito_candidate.fasta
```
then used blast and mitos against my sequence and found duplicate hits across the 23kb region 
using positions from mitos output i can remove the duplicates of the mitochondria (duplicates occur because its circular) 

```
seqkit subseq whitaker_mito_candidate.fasta -r 1428:18720 -o whitaker_mito_trimmed_full.fasta
```



blasting mitos annotations 
convert bed to fa 
```
bedtools getfasta \
    -fi whitaker_mito_trimmed_full.fasta \
    -bed mitos_annotations.bed \
    -s \
    -name \
    -fo mitos_genes.fa
```

then blast individual genes 
```
blastn -query mitos_genes.fa \
       -db nt \
       -remote \
       -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
       -max_target_seqs 5 \
       -out genes_vs_nt.tsv
```






















# Generating my genbank file

whitaker_mito_to_gb.py
```
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Input files
fasta_file = "whitaker_mito_trimmed_full.fasta"
bed_file = "mitos_annotations.bed"
output_file = "whitaker_mito.gb"

# Load the FASTA sequence
record = SeqIO.read(fasta_file, "fasta")

# Create a list to hold features
features = []

# Parse BED
with open(bed_file) as bed:
    for line in bed:
        if line.startswith("#") or line.strip() == "":
            continue
        parts = line.strip().split("\t")
        start = int(parts[1])  # BED start (0-based)
        end = int(parts[2])    # BED end
        gene_name = parts[3] if len(parts) > 3 else "unknown"
        strand = 1
        if len(parts) > 5 and parts[5] == "-":
            strand = -1
        
        # Add feature
        feature = SeqFeature.SeqFeature(
            FeatureLocation(start, end, strand=strand),
            type="gene",
            qualifiers={"gene": gene_name}
        )
        features.append(feature)

# Add features to the record
record.features = features

# Write to GenBank
SeqIO.write(record, output_file, "genbank")

print(f"GenBank file written to {output_file}")
```


# Coverage 

align Q20 reads back to mitogenome fasta file to produce summary statistics of read depth
```
#!/bin/bash -e
#SBATCH --job-name=mito_cov
#SBATCH --account=uoo04250
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=mito_cov_%j.out
#SBATCH --error=mito_cov_%j.err

module --force purge
module load NeSI
module load minimap2/2.28-GCC-12.3.0
module load SAMtools/1.22-GCC-12.3.0

REF="whitaker_mitogenomenew.fasta"
READS="reads_q20_l1kv2.fastq.gz"
PREFIX="whitaker_mitonewcoverage"

minimap2 -ax map-ont -t ${SLURM_CPUS_PER_TASK} "$REF" "$READS" | \
    samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${PREFIX}.sorted.bam

samtools index ${PREFIX}.sorted.bam
samtools depth -a ${PREFIX}.sorted.bam > ${PREFIX}.depth.txt

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
```

Mean_coverage	527.1
Min_coverage	79
Max_coverage	693
Pct_bases_>=1x	100
Pct_bases_>=10x	100
Pct_bases_>=20x	100
