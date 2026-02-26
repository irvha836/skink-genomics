Calculating genome wide heterozygosity with masked genome aligning q20 reads with a minimum length of 1000 bp

mpileup step to generate variant sites (all-sites VCF)

```
#!/bin/bash
#SBATCH --job-name=ONT_mpileup
#SBATCH --account=uoo04250
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=ONT_mpileup_%j.out
#SBATCH --error=ONT_mpileup_%j.err

# Load modules
module load SAMtools
module load BCFtools

# Paths
REF="/nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb/whitakergenome_hardmasked.fasta"
BAM="/nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb/reads.sorted.q20.bam"
OUT_DIR="/nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb"
VCF="$OUT_DIR/reads.allsites.vcf.gz"

# Make output directory
mkdir -p $OUT_DIR

# Step 1: Index reference if not done already
if [ ! -f "$REF.fai" ]; then
    echo "Indexing reference..."
    samtools faidx $REF
fi

# Step 2: Run bcftools mpileup → call all sites
echo "Generating all-sites VCF..."
bcftools mpileup \
    -f $REF \
    -a AD,DP \
    -Q 10 -q 20 \
    -Ou \
    $BAM \
| bcftools call -m -Oz -o $VCF

# Step 3: Index VCF
echo "Indexing VCF..."
bcftools index $VCF

echo "All-sites VCF complete: $VCF"
```

filter vcf for low quality sites, but also keeps monomorphic sites 

investgiate depths filters

```
module load mosdepth
mosdepth depth reads.sorted.q20.bam
```

```
awk 'NR>1 && $2>=10000 && $6<1000 {sum += $2*$4; len += $2} END {print "Genome-wide mean depth =", sum/len}' depth.mosdepth.summary.txt 
```

Genome-wide mean depth = 16.757
minimum depth 1/3rd mean depth =6
maximum depth 2x mean depth =34

```
module load BCFtools

bcftools filter \
-e 'DP<6 || DP>34' \
/nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb/reads.allsites.vcf.gz \
-Oz \
-o /nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb/reads.allsites.DPfiltered.vcf.gz

bcftools index /nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb/reads.allsites.DPfiltered.vcf.gz

```

generate contig lengths for per window/contig heterozygosity caluclations 
```
cut -f1,2 \
/nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb/whitakergenome_hardmasked.fasta.fai \
> /nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb/contig_lengths.txt

```




# **Streamlined Pipeline**

index and align q20 reads with minimum 1k bp in length to hardmasekd genome
```
#!/bin/bash
#SBATCH --account=uoo04250
#SBATCH --job-name=ONT_align
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=ONT_align_%j.out
#SBATCH --error=ONT_align_%j.err

# Load modules
module load minimap2
module load SAMtools

# Paths
REF="/nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb/whitakergenome_hardmasked.fasta"
READS="/nesi/nobackup/uoo04250/genome_assembly/psmcclean/reads_q20_l1kv2.fastq.gz"
OUT_DIR="/nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb"
BAM="$OUT_DIR/reads.sorted.q20.bam"

# Make output directory if it doesn't exist
mkdir -p $OUT_DIR

# Step 1: Build minimap2 index if not already done
if [ ! -f "$REF.mmi" ]; then
    echo "Building minimap2 index for reference..."
    minimap2 -d $REF.mmi $REF
fi

# Step 2: Align ONT reads → sort → filter low-quality reads
echo "Aligning reads and filtering..."
minimap2 -ax map-ont -t 8 $REF.mmi $READS \
    | samtools sort -@ 8 \
    | samtools view -@ 8 -b -q 20 > $BAM

# Step 3: Index the final BAM
echo "Indexing BAM..."
samtools index $BAM

echo "Alignment complete. BAM located at $BAM"
```

Generate all sites vcf with sorted reads 
```
module load SAMtools/1.16 BCFtools/1.16

bcftools mpileup \
    -f whitakergenome_hardmasked.fasta \
    -a AD,DP \
    -Q 10 -q 20 \
    -Ou reads.sorted.q20.bam \
| bcftools call -m --threads 8 -Oz -o reads.allsites.vcf.gz

```
investigate depth filter 

```
module load mosdepth
mosdepth depth reads.sorted.q20.bam
```
```
awk 'NR>1 && $2>=10000 && $6<1000 {sum += $2*$4; len += $2} END {print "Genome-wide mean depth =", sum/len}' depth.mosdepth.summary.txt 
```
Genome-wide mean depth = 16.757
minimum depth 1/3rd mean depth =6
maximum depth 2x mean depth =34

filter for depth and generate chromosome lengths
```
#!/bin/bash
#SBATCH --job-name=ONT_filter_only
#SBATCH --account=uoo04250
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=ONT_filter_only_%j.out
#SBATCH --error=ONT_filter_only_%j.err

module purge
module load SAMtools/1.16
module load BCFtools/1.16

OUT_DIR="/nesi/nobackup/uoo04250/genome_assembly/genome_het_1stfeb"
REF="$OUT_DIR/whitakergenome_hardmasked.fasta"
VCF_ALL="$OUT_DIR/reads.allsites.vcf.gz"
VCF_FILTERED="$OUT_DIR/reads.allsites.DPfiltered.vcf.gz"
CHROM_LEN="$OUT_DIR/chromosome_lengths.txt"

# Step 3: Apply depth filter (FIXED)
bcftools filter \
    -e 'INFO/DP<6 || INFO/DP>34' \
    -Oz \
    -o $VCF_FILTERED \
    $VCF_ALL

# Step 4: Index filtered VCF
bcftools index $VCF_FILTERED

# Step 5: Generate contig lengths from FASTA index
cut -f1,2 ${REF}.fai > $OUT_DIR/contig_lengths_raw.txt
awk '$2>=10000 {print $1 "\t" $2}' \
    $OUT_DIR/contig_lengths_raw.txt > $CHROM_LEN

echo "Filtering and contig length generation complete."
```

generate genome-wide heterozygosity by calculating it per contig with python script
```
#!/bin/bash -e
#SBATCH --job-name=gwh_allcontigs
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --account=uoo04250
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --hint=nomultithread

module purge
module load Python/3.11.6-foss-2023a

VCF_FILE="reads.allsites.DPfiltered.vcf.gz"
CHROM_LENGTHS="chromosome_lengths.txt"
CONTIGS="contigs.txt"
WINDOW_SIZE=1000000
STEP_SIZE=100000

while read -r chrom; do
    echo "Processing contig: $chrom"
    python ./Heterozygosity_windows.py \
        "$VCF_FILE" \
        "$CHROM_LENGTHS" \
        "$WINDOW_SIZE" \
        "$STEP_SIZE" \
        "$chrom"
done < "$CONTIGS"

```

ONT reads → single-sample genotyping → depth filtering critical

Hard-masked reference reduces false positives in repetitive regions

PASS filtering ensures only high-confidence sites are used

Contigs <10 kb excluded to avoid tiny, noisy windows

Window size / step size can be adjusted depending on genome size and resolution needed

code adapted from sebastian and from  https://www.science.org/doi/full/10.1126/sciadv.aau0757 
