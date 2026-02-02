use transition/transversion ratio of 2.0 from https://academic.oup.com/g3journal/article/12/2/jkab402/6433156 if bam fails. 



Calculating genome wide heterozygosity with masked genome aligning q20 reads with a minimum length of 1000 bp

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

filter vcf for low quality sites but also keeps monomorphic sites 

investgiate depths filters

```
module load mosdepth
mosdepth depth reads.sorted.q20.bam

```
awk 'NR>1 && $2>=10000 && $6<1000 {sum += $2*$4; len += $2} END {print "Genome-wide mean depth =", sum/len}' depth.mosdepth.summary.txt 

```

Genome-wide mean depth = 16.757
module load BCFtools

minimum depth 1/3rd mean depth =6
maximum depth 2x mean depth =34
```
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
