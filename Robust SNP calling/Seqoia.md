Generate a new popmap with harsh filtering (r80 and max-obs-het-65%)for better estimates, aiming for 300-700 SNPs my castle popmap excluded low coverage individuals C05, C06, C07 and included the legacy sample VUW05

# Castle individuals
```
mkdir -p castle_population_filtered
populations -P gstacks_out/ \
-M castle_popmap.txt \
-O castle_population_filtered/ \
--vcf \
--max-obs-het 0.65 \
-R 0.8
```
```
module load VCFtools
vcftools --vcf populations_out/populations.snps.vcf --missing-indv
sort -k 5n out.imiss > out.imiss.sorted
```

Further filter for minDP of 5 and max missing of 0.8 
```
vcftools --vcf populations.snps.vcf \
  --minDP 5 \
  --max-missing 0.8 \
  --recode \
  --out castle_filtered_minDP5_maxmiss80
```
gave me 169 snps out of 4110 sites

lowering mindp to 3
```
vcftools --vcf populations.snps.vcf \
  --minDP 3 \
  --max-missing 0.8 \
  --recode \
  --out castle_filtered_minDP3_maxmiss80
```
gave me 1111 out of a possible 4110 Sites

mindp4 kept 474 out of a possible 4110 Sites

to get snp counts 
```
bcftools view -H castle_filtered_minDP3_maxmiss80.recode.vcf | wc -l
```

# Northland individuals

removed castle and admixed individuals from popmap as well as individuals with more than 90% missing data apart from pedigree important individuals such as F21 & F23 
(D01, F09, B01, A02, D09, D29, F02, F28, A03, F13, F10, A04, D04, F25, DN75, PS44)

```
populations -P gstacks_out/ \
-M northland_popmap.txt \
-O northland_population_filtered/ \
--vcf \
--max-obs-het 0.65 \
-R 0.8
```
```
vcftools --vcf populations.snps.vcf \
  --minDP 5 \
  --max-missing 0.8 \
  --recode \
  --out northland_filtered_minDP5_maxmiss80
 ```

Keeping 57 out of 57 individuals. However, this only retained 5 out of 6898 snps.


Tried again with Mindp3 max-missing 0.8 kept 273 out of a possible 6898 Sites, mindp 2 kept 1283 out of a possible 6898 Sites.

#PLINK
module load PLINK/1.09b6.16

plink \
  --vcf northland_filtered_minDP2_maxmiss80.recode.vcf \
  --allow-extra-chr \
  --recode A \
  --out northland















# New Workflow after talking to Yasmin 
she recomended no LD purging or maf however since my depth is low and SNPs are not great i decided to go ahead anyway with relaxed paramaters to try get more reliable SNPs to infer the pedigree. 
```
vcftools --vcf populations.snps.vcf \
  --minDP 2 \
  --maf 0.01 \
  --max-missing 0.8 \
  --recode \
  --out northland_filtered_maf001_mm80_mindp2
```
After filtering, kept 1189 out of a possible 6898 Sites
Run Time = 0.00 seconds

```
module load PLINK/1.09b6.16

plink --vcf northland_filtered_maf001_mm80_mindp2.recode.vcf  \
  --allow-extra-chr \
  --make-bed \
  --out northlandv2
```
prune for linkage disequilibruim relaxed (r2=0.3) 
```
plink --bfile northland \
      --allow-extra-chr \
      --indep-pairwise 50 5 0.3 \
      --out northland_pruned
```
885 SNPs left after LD pruning

extract pruned SNPs for seqoia 
```
plink --bfile northlandv2 \
      --allow-extra-chr \
      --extract northland_prunedv2.prune.in \
      --recode A \
      --out northland_sequoia_v2
```
download raw output file and use for sequoia 

# NO LD 
```
plink --vcf northland_filtered_maf001_mm80_mindp2.recode.vcf  \
  --allow-extra-chr \
  --make-bed \
  --out northland_no_LD
```
```
plink --bfile northland_no_LD \
 --allow-extra-chr \
 --recode A  \
 --out northland_no_LD_sequoia
```
retained 1189 variant sites.


# Sequoia R
setwd("~/Downloads")

install.packages("sequoia")
library(sequoia) 
library(data.table)

plink_df <- read.table("northland_sequoia_v2.raw", header=TRUE)
head(plink_df)
colnames(plink_df)[1:10]

GenoM <- GenoConvert(InFile="northland_sequoia_v2.raw",
                     InFormat="raw")


#Read in the LifeHistData file
LH <- read.table("LifeHistData.txt", 
                 header = TRUE,   # your file has column names
                 sep = "\t",      # tab-delimited
                 stringsAsFactors = FALSE)  # keep character columns as characters


#Read in the PLINK raw file
plink_df <- read.table("northland_sequoia_v2.raw", header = TRUE)

#Check first few IDs
head(plink_df$IID)

#Replace IDs starting with "F" or "M" that clash with Sequoia dummy substitutions
#Example: add "H" prefix instead
#(I adjust the IDs of F01-F27)
plink_df$IID <- gsub("^F", "H", plink_df$IID)
plink_df$IID <- gsub("^M", "H", plink_df$IID)

#Check the new IDs
head(plink_df$IID)

#Save to a new safe raw file
write.table(plink_df, "northland_sequoia_v2_safe.raw",
            quote = FALSE, row.names = FALSE, sep = "\t")

#Now run GenoConvert safely
GenoM <- GenoConvert(InFile = "northland_sequoia_v2_safe.raw",
                     InFormat = "raw")


plink_df <- read.table("northland_sequoia_v2_safe.raw", header = TRUE)
plink_df_dt <- as.data.table(plink_df)
cols_to_convert <- grep("^X", names(plink_df_dt), value = TRUE)  # X = your SNP columns
plink_df_dt[, (cols_to_convert) := lapply(.SD, as.character), .SDcols = cols_to_convert]

#convert genotypes to seqouia format
GenoM <- GenoConvert(InFile = plink_df_dt, InFormat = 'raw')
GenoM.checked <- CheckGeno(GenoM, Return = "GenoM")

LH <- read.table("LifeHistData.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ParOUT <- sequoia(GenoM = GenoM, LifeHistData = LH,
                  Module="par", Err=0.01, quiet=FALSE, Plot=TRUE)
ParOUT$DupGenotype  # check for duplicates

IndivMiss <- apply(GenoM, 1, function(x) sum(x == -9) / ncol(GenoM))
IndivMiss[IndivMiss > 0.8]  # anyone with <20% SNPs


SeqOUT <- sequoia(GenoM = GenoM,
                  LifeHistData = LH,
                  Module = "ped",   # full pedigree
                  Err = 0.01)       # SNP genotyping error rate

#Summary of pedigree
SummarySeq(SeqOUT)

#Check the assigned parents for H21
SeqOUT$PedigreePar[SeqOUT$PedigreePar$id == "H21", ]

#Check full pedigree (includes dummy parents)
SeqOUT$Pedigree[SeqOUT$Pedigree$id == "H21", ]

write.csv(SeqOUT$PedigreePar, file = "PedigreeParentage.csv", row.names = FALSE)
write.csv(SeqOUT$Pedigree, file = "Pedigree_full.csv", row.names = FALSE)
write.csv(SeqOUT$DummyIDs, file = "Pedigree_dummyIDs.csv", row.names = FALSE)


stats <- SnpStats(GenoM, SeqOUT$PedigreePar)


Maybe <- GetMaybeRel(
  GenoM = GenoM,
  Pedigree = SeqOUT$PedigreePar,
  LifeHistData = LH,
  Err = 0.01,
  Complex = "full",
  Module = "ped"
)

#Inspect
Maybe$MaybeRel

#Save
write.csv(Maybe$MaybeRel, file = "MaybeRelatives.csv", row.names = FALSE)


MaybeM <- GetRelM(Pairs = Maybe$MaybeRel)
PlotRelPairs(MaybeM)


RelM <- GetRelM(Pedigree = SeqOUT$PedigreePar, Pairs = Maybe$MaybeRel)
PlotRelPairs(RelM)


#Pedigree plot including parental assignment
Rel.sd <- GetRelM(SeqOUT$PedigreePar, patmat = TRUE, GenBack = 2)
PlotRelPairs(Rel.sd)

#Full pedigree including dummy parents
Rel.sd2 <- GetRelM(SeqOUT$Pedigree, patmat = TRUE, GenBack = 2)
PlotRelPairs(Rel.sd2)

