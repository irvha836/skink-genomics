
#plotting relatedness
Sys.which("gfortran")
system("R CMD config FC", intern = TRUE)

install.packages("related", repos="https://R-Forge.R-project.org", type="source")
library(related)
setwd("~/Downloads")
library(related)
library(data.table)

# read plink raw file
geno <- fread("northland_no_LD_sequoia.raw")

# keep IDs
ids <- geno$IID

# extract SNP dosage columns only
X <- as.data.frame(geno[, 7:ncol(geno), with = FALSE])

# convert 0/1/2 plink dosage into two allele columns per SNP
split_snp <- function(x) {
  a1 <- ifelse(x == 2, 2, 1)
  a2 <- ifelse(x == 0, 1, 2)
  a1[is.na(x)] <- 0
  a2[is.na(x)] <- 0
  out <- data.frame(a1, a2)
  return(out)
}

geno2 <- do.call(cbind, lapply(X, split_snp))

# add IDs as first column
geno_related <- data.frame(ID = ids, geno2, check.names = FALSE)

# run TrioML
res <- coancestry(
  geno = geno_related,
  trioml = 1,
  wang = 0,
  lynchli = 0,
  lynchrd = 0,
  quellergt = 0,
  ritland = 0,
  dyadml = 0,
  allow.inbreeding = TRUE
)

# view results
head(res$relatedness)

# save results
write.csv(res$relatedness, "northland_trioml_results.csv", row.names = FALSE)

#Plotting
library(dplyr)
library(ggplot2)
library(reshape2)

rel <- read.csv("northland_trioml_results.csv", stringsAsFactors = FALSE)

ids <- sort(unique(c(rel$ind1.id, rel$ind2.id)))

rel_mat <- matrix(
  NA_real_,
  nrow = length(ids),
  ncol = length(ids),
  dimnames = list(ids, ids)
)

for (i in seq_len(nrow(rel))) {
  id1 <- rel$ind1.id[i]
  id2 <- rel$ind2.id[i]
  val <- rel$trioml[i]
  
  rel_mat[id1, id2] <- val
  rel_mat[id2, id1] <- val
}

# keep self-relatedness on diagonal
diag(rel_mat) <- 0.5

# replace any remaining missing values with 0
rel_mat[is.na(rel_mat)] <- 0

# cluster individuals
hc <- hclust(dist(rel_mat))
rel_mat <- rel_mat[hc$order, hc$order]

# remove only lower triangle, keep diagonal
rel_mat[lower.tri(rel_mat, diag = FALSE)] <- NA

rel_long <- melt(rel_mat, varnames = c("Ind1", "Ind2"), value.name = "trioml")

p <- ggplot(rel_long, aes(x = Ind2, y = Ind1, fill = trioml)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradientn(
    colours = c("#f6f6d5","#c2e699","#78c679","#2c7fb8","#253494"),
    limits = c(0,0.7),
    breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),
    na.value = "white",
    name = expression("Kinship coefficient (" * phi * ")")
  ) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    legend.key.height = unit(2, "cm"),
    legend.key.width  = unit(0.6, "cm")
  )

summary(rel$trioml)
max(rel$trioml, na.rm = TRUE)
rel %>% filter(trioml > 0.5)

pdf("TRIOML.pdf", width = 6, height = 5, family = "Times", useDingbats = FALSE)
print(p)
dev.off()

pdf("TRIOML.pdf", width = 6, height = 5, family = "Times", useDingbats = FALSE)
print(p)
dev.off()







#plotting kinship r/2
Sys.which("gfortran")
system("R CMD config FC", intern = TRUE)

install.packages("related", repos = "https://R-Forge.R-project.org", type = "source")
library(related)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)

setwd("~/Downloads")

# read plink raw file
geno <- fread("northland_no_LD_sequoia.raw")

# keep IDs
ids <- geno$IID

# extract SNP dosage columns only
X <- as.data.frame(geno[, 7:ncol(geno), with = FALSE])

# convert 0/1/2 plink dosage into two allele columns per SNP
split_snp <- function(x) {
  a1 <- ifelse(x == 2, 2, 1)
  a2 <- ifelse(x == 0, 1, 2)
  a1[is.na(x)] <- 0
  a2[is.na(x)] <- 0
  data.frame(a1, a2)
}

geno2 <- do.call(cbind, lapply(X, split_snp))

# add IDs as first column
geno_related <- data.frame(ID = ids, geno2, check.names = FALSE)

# run TrioML
res <- coancestry(
  geno = geno_related,
  trioml = 1,
  wang = 0,
  lynchli = 0,
  lynchrd = 0,
  quellergt = 0,
  ritland = 0,
  dyadml = 0,
  allow.inbreeding = TRUE
)

# save results
write.csv(res$relatedness, "northland_trioml_results.csv", row.names = FALSE)

# ----------------------------
# Plot TrioML as KINSHIP (phi)
# ----------------------------
rel <- read.csv("northland_trioml_results.csv", stringsAsFactors = FALSE)

ids <- sort(unique(c(rel$ind1.id, rel$ind2.id)))

rel_mat <- matrix(
  NA_real_,
  nrow = length(ids),
  ncol = length(ids),
  dimnames = list(ids, ids)
)

for (i in seq_len(nrow(rel))) {
  id1 <- rel$ind1.id[i]
  id2 <- rel$ind2.id[i]
  val_r <- rel$trioml[i]       # TrioML relatedness
  val_phi <- val_r / 2         # convert relatedness to kinship
  
  rel_mat[id1, id2] <- val_phi
  rel_mat[id2, id1] <- val_phi
}

# self-kinship on phi scale
diag(rel_mat) <- 0.5

# replace remaining missing values with 0
rel_mat[is.na(rel_mat)] <- 0

# cluster individuals using phi matrix
hc <- hclust(dist(rel_mat))
rel_mat <- rel_mat[hc$order, hc$order]

# keep upper triangle
rel_mat[lower.tri(rel_mat, diag = FALSE)] <- NA

rel_long <- melt(rel_mat, varnames = c("Ind1", "Ind2"), value.name = "phi")

p <- ggplot(rel_long, aes(x = Ind2, y = Ind1, fill = phi)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradientn(
    colours = c("#f6f6d5", "#c2e699", "#78c679", "#2c7fb8", "#253494"),
    limits = c(0, 0.5),
    breaks = c(0, 0.0625, 0.125, 0.25, 0.5),
    na.value = "white",
    name = expression("Kinship coefficient (" * phi * ")")
  ) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    legend.key.height = unit(2, "cm"),
    legend.key.width  = unit(0.6, "cm")
  )

print(p)

pdf("TRIOML_kinship_phi.pdf", width = 6, height = 5, family = "Times", useDingbats = FALSE)
print(p)
dev.off()

ggsave(
  "TRIOMLkinships.pdf",
  plot = p,
  width = 6,
  height = 5,
  device = cairo_pdf,
  family = "Times New Roman"
)

pdf("TRIOMLkinships.pdf", width = 6, height = 5, family = "Times", useDingbats = FALSE)
print(p)
dev.off()
# ----------------------------
# Useful summaries
# ----------------------------
summary(rel$trioml)                  # original TrioML relatedness scale
max(rel$trioml, na.rm = TRUE)

# highest raw TrioML relatedness pairs
rel %>%
  filter(ind1.id != ind2.id) %>%
  mutate(phi_trioml = trioml / 2) %>%
  arrange(desc(trioml)) %>%
  select(ind1.id, ind2.id, trioml, phi_trioml) %>%
  slice_head(n = 800)



#correlation between seq and trio coefficent
library(data.table)

# ----------------------------
# 1. Rebuild TrioML matrix fresh as KINSHIP (phi)
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
  val <- rel_df$trioml[i] / 2   # convert TrioML relatedness r to kinship phi
  
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
cor.test(df$sequoia, df$trioml, method = "pearson")

# only Sequoia-inferred relationships
df_seq_rel <- df[df$sequoia > 0, ]
cor.test(df_seq_rel$sequoia, df_seq_rel$trioml, method = "spearman")



df$seq_related <- df$sequoia > 0
df$tri_related <- df$trioml > 0.025

tab <- table(df$seq_related, df$tri_related)
tab
