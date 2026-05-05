setwd("~/Downloads")
library(sequoia)
library(data.table)
library(dplyr)

# --- Read LifeHistData ONCE and apply same rename as genotypes (F -> H)
LH <- read.table("LifeHistData.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
LH$ID <- sub("^F", "H", LH$ID)

# --- Read PLINK raw, rename IID F->H, write safe file
plink_df <- read.table("northland_no_LD_sequoia.raw", header=TRUE)
plink_df$IID <- sub("^F", "H", plink_df$IID)   # (keep the ^M line removed since you don't need it)

write.table(plink_df, "northland_no_LD_sequoia_safe.raw",
            quote=FALSE, row.names=FALSE, sep="\t")

# --- Convert genotypes
plink_df_dt <- as.data.table(read.table("northland_no_LD_sequoia_safe.raw", header=TRUE))
cols_to_convert <- grep("^X", names(plink_df_dt), value=TRUE)
plink_df_dt[, (cols_to_convert) := lapply(.SD, as.character), .SDcols = cols_to_convert]

GenoM <- GenoConvert(InFile = plink_df_dt, InFormat = "raw")
GenoM <- CheckGeno(GenoM, Return="GenoM")

# --- Run sequoia
ParOUT <- sequoia(GenoM=GenoM, LifeHistData=LH, Module="par", Err=0.01, Plot=TRUE)

SeqOUT <- sequoia(GenoM=GenoM, LifeHistData=LH, Module="ped", Err=0.01)
#up to here
SummarySeq(SeqOUT)

SeqOUT$PedigreePar[SeqOUT$PedigreePar$id == "H21", ]

# Check full pedigree (includes dummy parents)
SeqOUT$Pedigree[SeqOUT$Pedigree$id == "H21", ]


Maybe <- GetMaybeRel(
  GenoM = GenoM,
  Pedigree = SeqOUT$PedigreePar,
  LifeHistData = LH,
  Err = 0.01,
  Complex = "full",
  Module = "ped"
)

# Inspect
Maybe$MaybeRel


MaybeM <- GetRelM(Pairs = Maybe$MaybeRel)
PlotRelPairs(MaybeM)


RelM <- GetRelM(Pedigree = SeqOUT$PedigreePar, Pairs = Maybe$MaybeRel)
PlotRelPairs(RelM)


# Pedigree plot including parental assignment
Rel.sd <- GetRelM(SeqOUT$PedigreePar, patmat = TRUE, GenBack = 2)
PlotRelPairs(Rel.sd)

pdf("RelPairs_plot.pdf",
    width = 7, height = 6,
    family = "Times")  # optional

Rel.sd <- GetRelM(SeqOUT$PedigreePar, patmat = TRUE, GenBack = 2)
PlotRelPairs(Rel.sd)

dev.off()

# Full pedigree including dummy parents
Rel.sd2 <- GetRelM(SeqOUT$Pedigree, patmat = TRUE, GenBack = 2)
PlotRelPairs(Rel.sd2)





#SAVE FILES
write.csv(SeqOUT$PedigreePar, "northlandnoLD_PedigreeParentage.csv", row.names = FALSE)
write.csv(SeqOUT$Pedigree,    "northlandnoLD_PedigreeFull.csv",      row.names = FALSE)
write.csv(SeqOUT$DummyIDs,    "northlandnoLD_DummyIDs.csv",          row.names = FALSE)

write.csv(Maybe$MaybeRel,     "northlandnoLD_MaybeRelatives.csv",    row.names = FALSE)


#QC
stats <- SnpStats(GenoM, SeqOUT$PedigreePar)
MAF <- ifelse(stats[,"AF"] <= 0.5, stats[,"AF"], 1 - stats[,"AF"])

summary(MAF)
summary(stats[,"Err.hat"])


with(SeqOUT$PedigreePar, table(sire == "0" | is.na(sire), dam == "0" | is.na(dam)))
nrow(SeqOUT$DummyIDs)
nrow(Maybe$MaybeRel)
table(Maybe$MaybeRel$TopRel)








colnames(LH)
table(LH$Sex, useNA="ifany")
summary(LH$BirthYear)
colnames(stats)








class(Rel.sd)
dim(Rel.sd)
Rel.sd[1:10, 1:10]
unique(as.vector(Rel.sd))
table(as.vector(Rel.sd), useNA = "ifany")







table(cut(phi, breaks = c(-0.01, 0.01, 0.1, 0.2, 0.3, 1)))




phi_no_diag <- phi
diag(phi_no_diag) <- NA

table(cut(phi_no_diag, breaks = c(-0.01, 0.01, 0.1, 0.2, 0.3, 1)))






# convert to kinship 

library(reshape2)
library(dplyr)
library(ggplot2)

rel <- as.matrix(Rel.sd)
rel <- apply(rel, c(1,2), as.character)

phi <- matrix(0, nrow = nrow(rel), ncol = ncol(rel),
              dimnames = dimnames(rel))

phi[rel == "S"] <- 0.5
phi[rel %in% c("P","O")] <- 0.25
phi[rel == "FS"] <- 0.25
phi[rel %in% c("MHS","PHS")] <- 0.125
phi[rel %in% c("PGF","PGM","MGF","MGM","GO")] <- 0.125
phi[rel %in% c("FA","FN")] <- 0.125
phi[rel %in% c("HA","HN")] <- 0.0625
phi[rel == "U"] <- 0

# make symmetric
phi <- pmax(phi, t(phi), na.rm = TRUE)

# TrioML order
trioml_order <- c(
  "D12","D10","D30","G03","G01","E05","D27","D19","D16","G02",
  "D28","D17","D15","E01","D14","D23","D18","D13","D11",
  "F19","F18","F14","D20","F21","F03","D05","F17","F11","F01",
  "F24","F16","D25","F07","D31","E03","D26","F06","F05","F23",
  "F22","F08","F04","D06","F15","F27","D22","E02","D21","D08",
  "D07","D03","D24","E04","F20","F12","D02","F26"
)

# convert F IDs to H IDs for sequoia
seqoia_order <- gsub("^F", "H", trioml_order)

# check before subsetting
setdiff(seqoia_order, rownames(phi))
setdiff(rownames(phi), seqoia_order)

# reorder phi to match TrioML order
phi <- phi[seqoia_order, seqoia_order]

# rename back to F for plotting
plot_names <- gsub("^H", "F", rownames(phi))
rownames(phi) <- plot_names
colnames(phi) <- plot_names

# 🔄 FLIP BOTH AXES (KEY STEP)
phi <- phi[rev(rownames(phi)), rev(colnames(phi))]

# keep mirrored half (lower-left visually)
phi[upper.tri(phi, diag = FALSE)] <- NA

# melt
phi_long <- melt(phi, varnames = c("Ind1", "Ind2"), value.name = "phi") %>%
  filter(!is.na(phi))

# set factor order (now already flipped)
phi_long$Ind1 <- factor(phi_long$Ind1, levels = rownames(phi))
phi_long$Ind2 <- factor(phi_long$Ind2, levels = colnames(phi))

# plot
p <- ggplot(phi_long, aes(x = Ind2, y = Ind1, fill = phi)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradientn(
    colours = c("#f6f6d5","#c2e699","#78c679","#2c7fb8","#253494"),
    limits = c(0, 0.5),
    breaks = c(0, 0.0625, 0.125, 0.25, 0.5),
    name = expression("Expected kinship (" * phi * ")")
  ) +
  coord_fixed() +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Times"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7)
  )

print(p)

pdf("seqoia.pdf", width = 6, height = 5, family = "Times", useDingbats = FALSE)
print(p)
dev.off()












#correlation between trio results and seqoia  
library(data.table)

# ----------------------------
# 1. Rebuild TrioML matrix fresh
# ----------------------------
rel_df <- fread("northland_trioml_results.csv")

ids <- sort(unique(c(rel_df$ind1.id, rel_df$ind2.id)))

trioml_full <- matrix(
  NA_real_,
  nrow = length(ids),
  ncol = length(ids),
  dimnames = list(ids, ids)
)

for (i in seq_len(nrow(rel_df))) {
  id1 <- rel_df$ind1.id[i]
  id2 <- rel_df$ind2.id[i]
  val <- rel_df$trioml[i]
  
  trioml_full[id1, id2] <- val
  trioml_full[id2, id1] <- val
}

diag(trioml_full) <- 0.5
trioml_full[is.na(trioml_full)] <- 0

# ----------------------------
# 2. Rebuild Sequoia phi matrix fresh
# ----------------------------
rel_seq <- as.matrix(Rel.sd)
rel_seq <- apply(rel_seq, c(1,2), as.character)

phi_full <- matrix(
  0,
  nrow = nrow(rel_seq),
  ncol = ncol(rel_seq),
  dimnames = dimnames(rel_seq)
)

phi_full[rel_seq == "S"] <- 0.5
phi_full[rel_seq %in% c("P","O")] <- 0.25
phi_full[rel_seq == "FS"] <- 0.25
phi_full[rel_seq %in% c("MHS","PHS")] <- 0.125
phi_full[rel_seq %in% c("PGF","PGM","MGF","MGM","GO")] <- 0.125
phi_full[rel_seq %in% c("FA","FN")] <- 0.125
phi_full[rel_seq %in% c("HA","HN")] <- 0.0625
phi_full[rel_seq == "U"] <- 0

phi_full <- pmax(phi_full, t(phi_full), na.rm = TRUE)

# rename H back to F so names match TrioML
rownames(phi_full) <- gsub("^H", "F", rownames(phi_full))
colnames(phi_full) <- gsub("^H", "F", colnames(phi_full))

# ----------------------------
# 3. Put both in the same order
# ----------------------------
trioml_order <- c(
  "D12","D10","D30","G03","G01","E05","D27","D19","D16","G02",
  "D28","D17","D15","E01","D14","D23","D18","D13","D11",
  "F19","F18","F14","D20","F21","F03","D05","F17","F11","F01",
  "F24","F16","D25","F07","D31","E03","D26","F06","F05","F23",
  "F22","F08","F04","D06","F15","F27","D22","E02","D21","D08",
  "D07","D03","D24","E04","F20","F12","D02","F26"
)

# keep only IDs present in both
common_ids <- trioml_order[
  trioml_order %in% rownames(trioml_full) &
    trioml_order %in% rownames(phi_full)
]

trioml_cmp <- trioml_full[common_ids, common_ids]
phi_cmp    <- phi_full[common_ids, common_ids]

# ----------------------------
# 4. Extract matching pairwise values
# ----------------------------
diag(trioml_cmp) <- NA
diag(phi_cmp) <- NA

trioml_vec <- trioml_cmp[upper.tri(trioml_cmp)]
phi_vec    <- phi_cmp[upper.tri(phi_cmp)]

df <- data.frame(
  trioml  = as.numeric(trioml_vec),
  sequoia = as.numeric(phi_vec)
)

df <- na.omit(df)

str(df)
nrow(df)

# ----------------------------
# 5. Correlation
# ----------------------------
cor.test(df$sequoia, df$trioml, method = "spearman")

# Optional Pearson too
cor.test(df$sequoia, df$trioml, method = "pearson")


#i think this one is important pairwise spearmans test
df_seq_rel <- df[df$sequoia > 0, ]

cor.test(df_seq_rel$sequoia, df_seq_rel$trioml, method = "spearman")

# ----------------------------
# 6. Agreement test: both above zero
# ----------------------------
df$seq_related <- df$sequoia > 0
df$tri_related <- df$trioml > 0.05   # adjust threshold if needed

tab <- table(df$seq_related, df$tri_related)
tab
fisher.test(tab)


ggplot(df_seq_rel, aes(x = factor(sequoia), y = trioml)) +
  geom_boxplot() +
  theme_classic()




library(vegan)
#check results for how report when also using spearman. 
mantel(
  as.dist(trioml_cmp),
  as.dist(phi_cmp),
  method = "spearman",
  permutations = 9999
)




sum(df$sequoia > 0)

sum(df$trioml > 0.05)







library(dplyr)

rel <- read.csv("northland_trioml_results.csv", stringsAsFactors = FALSE)

# remove self-comparisons just in case
rel_pairs <- rel %>%
  filter(ind1.id != ind2.id) %>%
  arrange(desc(trioml))

# top 20 highest pairwise estimates
top25 <- rel_pairs %>%
  select(ind1.id, ind2.id, trioml) %>%
  slice_head(n = 25)

top25





#correlation trial two 

library(data.table)
library(vegan)

# ----------------------------
# 1. Rebuild TrioML matrix (convert r → kinship φ by dividing by 2)
# ----------------------------
rel_df <- fread("northland_trioml_results.csv")

ids <- sort(unique(c(rel_df$ind1.id, rel_df$ind2.id)))

trioml_full <- matrix(
  NA_real_,
  nrow = length(ids),
  ncol = length(ids),
  dimnames = list(ids, ids)
)

for (i in seq_len(nrow(rel_df))) {
  id1 <- rel_df$ind1.id[i]
  id2 <- rel_df$ind2.id[i]
  
  # convert relatedness (r) to kinship (φ)
  val <- rel_df$trioml[i] / 2
  
  trioml_full[id1, id2] <- val
  trioml_full[id2, id1] <- val
}

diag(trioml_full) <- 0.5
trioml_full[is.na(trioml_full)] <- 0

# ----------------------------
# 2. Rebuild Sequoia kinship matrix
# ----------------------------
rel_seq <- as.matrix(Rel.sd)
rel_seq <- apply(rel_seq, c(1, 2), as.character)

phi_full <- matrix(
  0,
  nrow = nrow(rel_seq),
  ncol = ncol(rel_seq),
  dimnames = dimnames(rel_seq)
)

phi_full[rel_seq == "S"] <- 0.5
phi_full[rel_seq %in% c("P", "O")] <- 0.25
phi_full[rel_seq == "FS"] <- 0.25
phi_full[rel_seq %in% c("MHS", "PHS")] <- 0.125
phi_full[rel_seq %in% c("PGF", "PGM", "MGF", "MGM", "GO")] <- 0.125
phi_full[rel_seq %in% c("FA", "FN")] <- 0.125
phi_full[rel_seq %in% c("HA", "HN")] <- 0.0625
phi_full[rel_seq == "U"] <- 0

# Force symmetry
phi_full <- pmax(phi_full, t(phi_full), na.rm = TRUE)

# Rename H → F so names match TrioML
rownames(phi_full) <- gsub("^H", "F", rownames(phi_full))
colnames(phi_full) <- gsub("^H", "F", colnames(phi_full))

# ----------------------------
# 3. Match matrices to same individuals + order
# ----------------------------
trioml_order <- c(
  "D12","D10","D30","G03","G01","E05","D27","D19","D16","G02",
  "D28","D17","D15","E01","D14","D23","D18","D13","D11",
  "F19","F18","F14","D20","F21","F03","D05","F17","F11","F01",
  "F24","F16","D25","F07","D31","E03","D26","F06","F05","F23",
  "F22","F08","F04","D06","F15","F27","D22","E02","D21","D08",
  "D07","D03","D24","E04","F20","F12","D02","F26"
)

common_ids <- trioml_order[
  trioml_order %in% rownames(trioml_full) &
    trioml_order %in% rownames(phi_full)
]

trioml_cmp <- trioml_full[common_ids, common_ids]
phi_cmp    <- phi_full[common_ids, common_ids]

# ----------------------------
# 4. Remove diagonal
# ----------------------------
diag(trioml_cmp) <- NA
diag(phi_cmp) <- NA

# ----------------------------
# 5. Mantel test — ALL pairs
# ----------------------------
mantel_all <- mantel(
  as.dist(trioml_cmp),
  as.dist(phi_cmp),
  method = "spearman",
  permutations = 9999,
  na.rm = TRUE
)

cat("\nMantel test: ALL pairs\n")
print(mantel_all)

# ----------------------------
# 6. Mantel test — RELATED pairs only
#    (TrioML φ > 0.025 OR Sequoia φ > 0)
# ----------------------------
threshold <- 0.025   # ⚠️ IMPORTANT: adjusted because we halved values

mask_related <- (trioml_cmp > threshold) | (phi_cmp > 0)

trioml_rel <- trioml_cmp
phi_rel    <- phi_cmp

trioml_rel[!mask_related] <- NA
phi_rel[!mask_related]    <- NA

mantel_related <- mantel(
  as.dist(trioml_rel),
  as.dist(phi_rel),
  method = "spearman",
  permutations = 9999,
  na.rm = TRUE
)

cat("\nMantel test: RELATED pairs only\n")
print(mantel_related)

# ----------------------------
# 7. Counts for reporting
# ----------------------------
upper_idx <- upper.tri(trioml_cmp)

n_all_pairs <- sum(!is.na(trioml_cmp[upper_idx]) & !is.na(phi_cmp[upper_idx]))
n_seq_related <- sum(phi_cmp[upper_idx] > 0, na.rm = TRUE)
n_tri_related <- sum(trioml_cmp[upper_idx] > threshold, na.rm = TRUE)
n_union_related <- sum(mask_related[upper_idx], na.rm = TRUE)

cat("\nCounts for reporting\n")
cat("All pairwise comparisons:", n_all_pairs, "\n")
cat("Sequoia-related pairs (phi > 0):", n_seq_related, "\n")
cat("TrioML-related pairs (> ", threshold, "): ", n_tri_related, "\n", sep = "")
cat("Related pairs (union):", n_union_related, "\n")

# ----------------------------
# 8. OPTIONAL: strict related pairs (both methods agree)
# ----------------------------
mask_strict <- (trioml_cmp > threshold) & (phi_cmp > 0)

trioml_strict <- trioml_cmp
phi_strict    <- phi_cmp

trioml_strict[!mask_strict] <- NA
phi_strict[!mask_strict]    <- NA

mantel_strict <- mantel(
  as.dist(trioml_strict),
  as.dist(phi_strict),
  method = "spearman",
  permutations = 9999,
  na.rm = TRUE
)

cat("\nMantel test: STRICT related pairs (both methods)\n")
print(mantel_strict)

n_strict_related <- sum(mask_strict[upper_idx], na.rm = TRUE)
cat("Strict related pairs:", n_strict_related, "\n")





# pairwise breeder kinship comparision between breeders 
library(dplyr)

# ----------------------------
# 1. Create breeder groups
# ----------------------------
ids <- rownames(trioml_cmp)

breeder <- case_when(
  grepl("^D", ids) ~ "D",
  grepl("^F", ids) ~ "F",
  grepl("^E", ids) ~ "E",
  grepl("^G", ids) ~ "G",
  TRUE ~ "Other"
)

names(breeder) <- ids

# ----------------------------
# 2. Extract pairwise values
# ----------------------------
upper_idx <- upper.tri(trioml_cmp)

tri_vals <- trioml_cmp[upper_idx]
seq_vals <- phi_cmp[upper_idx]

id1 <- rownames(trioml_cmp)[row(trioml_cmp)[upper_idx]]
id2 <- colnames(trioml_cmp)[col(trioml_cmp)[upper_idx]]

df_pairs <- data.frame(
  id1 = id1,
  id2 = id2,
  trioml = tri_vals,
  sequoia = seq_vals
)

# Add breeder info
df_pairs$group1 <- breeder[df_pairs$id1]
df_pairs$group2 <- breeder[df_pairs$id2]

# ----------------------------
# 3. Define within vs between
# ----------------------------
df_pairs <- df_pairs %>%
  mutate(
    within_group = group1 == group2
  )

# ----------------------------
# 4. Summary stats
# ----------------------------
df_pairs %>%
  group_by(within_group) %>%
  summarise(
    mean_trioml = mean(trioml, na.rm = TRUE),
    mean_sequoia = mean(sequoia, na.rm = TRUE),
    n = n()
  )


wilcox.test(trioml ~ within_group, data = df_pairs)
wilcox.test(sequoia ~ within_group, data = df_pairs)




#sequoia within and between 
df_seq <- df_pairs %>%
  filter(!is.na(sequoia))

df_seq %>%
  group_by(within_group) %>%
  summarise(
    mean_sequoia = mean(sequoia, na.rm = TRUE),
    n = n()
  )



df_seq_rel <- df_pairs %>%
  filter(sequoia > 0)

df_seq_rel %>%
  group_by(within_group) %>%
  summarise(
    mean_sequoia = mean(sequoia, na.rm = TRUE),
    n = n()
  )



# are trio estimates higher than sequoia 
wilcox.test(df_pairs$trioml, df_pairs$sequoia, paired = TRUE)

mask_strict <- (trioml_cmp > threshold) & (phi_cmp > 0)
tri_strict <- trioml_cmp
seq_strict <- phi_cmp

tri_strict[!mask_strict] <- NA
seq_strict[!mask_strict] <- NA
upper_idx <- upper.tri(tri_strict)

tri_vals <- tri_strict[upper_idx]
seq_vals <- seq_strict[upper_idx]
valid <- !is.na(tri_vals) & !is.na(seq_vals)

tri_vals <- tri_vals[valid]
seq_vals <- seq_vals[valid]
wilcox.test(tri_vals, seq_vals, paired = TRUE)

length(tri_vals)









wilcox.test(sequoia ~ within_group, data = df_seq_rel)



Rel.sd <- GetRelM(SeqOUT$PedigreePar, patmat = TRUE, GenBack = 2)

# rename individuals: H → F
new_names <- gsub("^H", "F", rownames(Rel.sd))

rownames(Rel.sd) <- new_names
colnames(Rel.sd) <- new_names

cairo_pdf("RelPairsF_plot.pdf",
          width = 7, height = 6)

PlotRelPairs(Rel.sd)

dev.off()
