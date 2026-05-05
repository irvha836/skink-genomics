library(ggplot2)
library(dplyr)

pca <- read.table("robust_PCA.eigenvec", header = FALSE)
colnames(pca)[1:2] <- c("FID","IID")
colnames(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))

eig <- scan("robust_PCA.eigenval")
var_exp <- eig / sum(eig) * 100


admixed_ids <- c(
  "A03","B01","A02","A01","G01","E05","G03",
  "D16","D27","D30","D19","D10","D12","D28","A04","B02"
)

# Create a 4-level grouping for plotting
pca$Group <- pca$Population
pca$Group[pca$IID %in% admixed_ids] <- "Admixed"

# (Optional) control legend order
pca$Group <- factor(pca$Group, levels = c("Castle","Northland","Wild","Admixed"))


#Mutora recode
library(ggplot2)
library(dplyr)

# Read PCA
pca <- read.table("robust_PCA.eigenvec", header = FALSE)
colnames(pca)[1:2] <- c("FID","IID")
colnames(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))

# Variance explained
eig <- scan("robust_PCA.eigenval")
var_exp <- eig / sum(eig) * 100

# Create grouping (Moturoa individual separated)
pca <- pca %>%
  mutate(Group = case_when(
    IID == "VUW06" ~ "Wild_Moturoa",
    grepl("^C", IID) ~ "Castle",
    grepl("^V", IID) ~ "Wild",
    TRUE ~ "Northland"
  ))

# Factor order for legend
pca$Group <- factor(pca$Group,
                    levels = c("Castle", "Northland", "Wild", "Wild_Moturoa"))

# Plot
p_PCA <- ggplot(pca, aes(PC1, PC2, fill = Group, shape = Group)) +
  geom_point(color = "black", stroke = 0.3, size = 3, alpha = 0.9) +
  scale_fill_manual(values = c(
    Castle = "red3",
    Northland = "#4575b4",
    Wild = "grey40",
    Wild_Moturoa = "black"   # distinct colour for VUW06
  )) +
  scale_shape_manual(values = c(
    Castle = 21,
    Northland = 21,
    Wild = 24,
    Wild_Moturoa = 24
  )) +
  theme_classic(base_size = 14, base_family = "Times") +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
    fill = "Population",
    shape = "Population"
  )

# Save
pdf("MutoroaPCA.pdf", width = 6, height = 5, family = "Times", useDingbats = FALSE)
print(p_PCA)
dev.off()
