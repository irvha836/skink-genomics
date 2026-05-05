exclude_ids <- c("B01", "B02", "A03", "A01", "A02", "A04", "VUW06")

# het
het <- read.table("robustpopulation.het", header=TRUE) %>%
  mutate(Ho = (N_SITES - O.HOM.) / N_SITES) %>%
  mutate(Population = case_when(
    grepl("^C", INDV) ~ "Castle",
    grepl("^V", INDV) ~ "Legacy",
    TRUE ~ "Northland"
  ))

# depth (change filename + column names as needed)
depth <- read.table("robustpopulation.idepth", header=TRUE)

# check names first:
colnames(depth)

# example rename (adjust to match your file)
depth <- depth %>%
  rename(mean_depth = MEAN_DEPTH)

# merge then exclude
dat <- het %>%
  left_join(depth, by = "INDV") %>%
  filter(!INDV %in% exclude_ids)

# now do the QC tests
summary(lm(Ho ~ mean_depth, data = dat))
summary(lm(Ho ~ Population + mean_depth, data = dat))

# and/or visualise
ggplot(dat, aes(mean_depth, Ho, color = Population)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic()

summary(lm(Ho ~ mean_depth, data = subset(dat, Population == "Northland")))
summary(lm(Ho ~ mean_depth, data = subset(dat, Population == "Legacy")))
summary(lm(Ho ~ mean_depth, data = subset(dat, Population == "Castle")))

head(data=depth INDV)

#Main Heterozygosity figure processing 
library(dplyr)
library(ggplot2)

admixed_ids <- c("A03","B01","A02","A01","G01","E05","G03",
                 "D16","D27","D30","D19","D10","D12","D28","A04","B02") 

# Read + compute Ho + assign populations + exclude admixed
het <- read.table("robustpopulation.het", header = TRUE) %>%
  mutate(
    Ho = (N_SITES - O.HOM.) / N_SITES,
    Population = case_when(
      grepl("^C", INDV) ~ "Castle",
      grepl("^V", INDV) ~ "Wild",
      TRUE ~ "Northland"
    )
  ) %>%
  filter(!INDV %in% exclude_ids) %>%
  mutate(
    Population = factor(Population, levels = c("Northland", "Legacy", "Castle"))
  )

# Heterozygosity summary table 
het %>%
  group_by(Population) %>%
  summarise(
    n = n(),
    mean_Ho = mean(Ho, na.rm = TRUE),
    sd_Ho = sd(Ho, na.rm = TRUE),
    median_Ho = median(Ho, na.rm = TRUE),
    min_Ho = min(Ho, na.rm = TRUE),
    max_Ho = max(Ho, na.rm = TRUE),
    .groups = "drop"
  )


#MAIN INBREEDING FIGURE 
library(dplyr)
library(ggplot2)

admixed_ids <- c("A03","B01","A02","A01","G01","E05","G03",
                 "D16","D27","D30","D19","D10","D12","D28")

fdat <- read.table("robustpopulation.het", header = TRUE) %>%
  mutate(
    Ho = (N_SITES - O.HOM.) / N_SITES,
    Group = case_when(
      INDV %in% admixed_ids ~ "Admixed",
      grepl("^C", INDV) ~ "Castle",
      grepl("^V", INDV) ~ "Wild",
      TRUE ~ "Northland"
    ),
    Group = factor(Group, levels = c("Wild", "Northland", "Castle", "Admixed"))
  )

# Summary table (inbreeding coefficient F)
f_summary <- fdat %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_F = mean(F, na.rm = TRUE),
    sd_F = sd(F, na.rm = TRUE),
    median_F = median(F, na.rm = TRUE),
    min_F = min(F, na.rm = TRUE),
    max_F = max(F, na.rm = TRUE),
    .groups = "drop"
  )

f_summary

p_F <- ggplot(fdat, aes(x = Group, y = F, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.70, linewidth = 0.8) +
  geom_jitter(aes(color = Group), width = 0.12, size = 1.8, alpha = 0.75) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2, fill = "white", stroke = 0.8) +
  scale_fill_manual(values = c(
    "Northland" = "#1f78b4",
    "Wild"    = "grey",
    "Castle"    = "#e31a1c",
    "Admixed"   = "purple"
  )) +
  scale_color_manual(values = c(
    "Northland" = "#1f78b4",
    "Wild"    = "grey30",
    "Castle"    = "#e31a1c",
    "Admixed"   = "purple4"
  )) +
  theme_classic(base_size = 15) + 
  labs(x = NULL, y = expression("Inbreeding coefficient (" * italic(F) * ")")) + 
  theme(
    legend.position = "none",
    text = element_text(family = "Times New Roman"),
    axis.line  = element_line(linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.9)
  )

p_F




fdat %>%
  filter(Group == "Admixed") %>%
  pull(INDV)

fdat %>%
  filter(Group == "Admixed") %>%
  arrange(F)

length(admixed_ids)

p_F <- ggplot(fdat, aes(x = Group, y = F, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.70, linewidth = 0.8) +
  geom_jitter(aes(color = Group), width = 0.12, size = 1.8, alpha = 0.75) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2, fill = "white", stroke = 0.8) +
  scale_fill_manual(values = c(
    "Northland" = "#1f78b4",
    "Wild" = "grey",
    "Castle" = "#e31a1c",
    "Admixed" = "purple"
  )) +
  scale_color_manual(values = c(
    "Northland" = "#1f78b4",
    "Wild" = "grey30",
    "Castle" = "#e31a1c",
    "Admixed" = "purple4"
  )) +
  theme_classic(base_size = 15) +
  labs(x = NULL, y = expression("Inbreeding coefficient (" * italic(F) * ")")) +
  theme(
    legend.position = "none",
    text = element_text(family = "Times"),
    axis.line  = element_line(linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.9)
  )



#save as pdf 
pdf("F_boxplotwildleftnoAo4.pdf", width = 6, height = 5, family = "Times", useDingbats = FALSE)
print(p_F)
dev.off()

#differences in inbreeding between the groups
fdat <- read.table("robustpopulation.het", header = TRUE) %>%
  mutate(
    Ho = (N_SITES - O.HOM.) / N_SITES,
    Group = case_when(
      INDV %in% admixed_ids ~ "Admixed",
      grepl("^C", INDV) ~ "Castle",
      grepl("^V", INDV) ~ "Wild",
      TRUE ~ "Northland"
    ),
    Group = factor(Group, levels = c("Northland", "Wild", "Castle", "Admixed"))
  )

kruskal.test(F ~ Group, data = fdat)

pairwise.wilcox.test(fdat$F, fdat$Group, p.adjust.method = "BH", exact = FALSE)




#Main heterozygosity figure 
library(dplyr)
library(ggplot2)

admixed_ids <- c("A03","B01","A02","A01","G01","E05","G03",
                 "D16","D27","D30","D19","D10","D12","D28",)

# Read + compute Ho + assign groups (including admixed)
het <- read.table("robustpopulation.het", header = TRUE) %>%
  mutate(
    Ho = (N_SITES - O.HOM.) / N_SITES,
    Group = case_when(
      INDV %in% admixed_ids ~ "Admixed",
      grepl("^C", INDV) ~ "Castle",
      grepl("^V", INDV) ~ "Wild",
      TRUE ~ "Northland"
    ),
    Group = factor(Group, levels = c("Wild", "Northland", "Castle", "Admixed"))
  )

# Heterozygosity summary table
het %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_Ho = mean(Ho, na.rm = TRUE),
    sd_Ho = sd(Ho, na.rm = TRUE),
    median_Ho = median(Ho, na.rm = TRUE),
    min_Ho = min(Ho, na.rm = TRUE),
    max_Ho = max(Ho, na.rm = TRUE),
    .groups = "drop"
  )

# Colours 
fill_cols <- c(
  "Northland" = "#1f78b4",
  "Wild"      = "grey",
  "Castle"    = "#e31a1c",
  "Admixed"   = "purple"
)

line_cols <- c(
  "Northland" = "#1f78b4",
  "Wild"      = "grey30",
  "Castle"    = "#e31a1c",
  "Admixed"   = "purple4"
)

# Main plot
p_Ho <- ggplot(het, aes(x = Group, y = Ho, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.70, linewidth = 0.8) +
  geom_jitter(aes(color = Group), width = 0.12, size = 1.8, alpha = 0.75) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2,
               fill = "white", stroke = 0.8) +
  scale_fill_manual(values = fill_cols) +
  scale_color_manual(values = line_cols) +
  theme_classic(base_size = 15, base_family = "Times") +
  labs(x = NULL, y = expression("Observed heterozygosity ("*H[o]*")")) +
  theme(
    legend.position = "none",
    axis.line  = element_line(linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.9)
  )

p_Ho

pdf("Ho_boxplotwildleftnoao4.pdf", width = 6, height = 5, family = "Times", useDingbats = FALSE)
print(p_Ho)
dev.off()


#testing correlation with percentage of admixture and heterozygosity 
q2 <- read.table("robust_captive_pruned.2.Q", header = FALSE)
colnames(q2) <- c("Cluster1","Cluster2")

fam <- read.table("robust_captive_pruned.fam")
colnames(fam) <- c("FID","IID","PAT","MAT","SEX","PHENO")

q2$IID <- fam$IID

admix <- q2 %>%
  rename(Castle = Cluster1, Northland = Cluster2)

het <- read.table("robustpopulation.het", header = TRUE) %>%
  mutate(Ho = (N_SITES - O.HOM.) / N_SITES)

het_admix <- het %>%
  left_join(admix, by = c("INDV" = "IID"))

#difference in Ho between groups
# Read + compute Ho + assign groups
het <- read.table("robustpopulation.het", header = TRUE) %>%
  mutate(
    Ho = (N_SITES - O.HOM.) / N_SITES,
    Group = case_when(
      INDV %in% admixed_ids ~ "Admixed",
      grepl("^C", INDV) ~ "Castle",
      grepl("^V", INDV) ~ "Wild",
      TRUE ~ "Northland"
    ),
    Group = factor(Group, levels = c("Northland", "Wild", "Castle", "Admixed"))
  )

# ---- Kruskal-Wallis test for heterozygosity ----
kruskal.test(Ho ~ Group, data = het)

pairwise.wilcox.test(het$Ho, het$Group, p.adjust.method = "BH", exact = FALSE)
#admixture correlation 
admixed <- het_admix %>%
  filter(Castle > 0.05 & Castle < 0.95)
cor.test(admixed$Ho, admixed$Castle)


p_admix <- ggplot(admixed, aes(x = Castle, y = Ho)) +
  geom_point(size = 3, color = "purple4", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  theme_classic(base_size = 15, base_family = "Times") +
  labs(
    x = "Castle ancestry proportion",
    y = expression("Observed heterozygosity ("*H[o]*")")
  )

p_admix

pdf("admixture_heterozygosity_correlation.pdf",
    width = 6, height = 5,
    family = "Times",
    useDingbats = FALSE)

print(p_admix)
dev.off()


#inbreeding correlation 

f_admix <- fdat %>%
  left_join(admix, by = c("INDV" = "IID"))

admixed_F <- f_admix %>%
  filter(Castle > 0.02 & Castle < 0.95)

cor.test(admixed_F$F, admixed_F$Castle)
cor.test(admixed_F$F, admixed_F$Castle, method = "pearson")

p_Fadmix <- ggplot(admixed_F, aes(x = Castle, y = F)) +
  geom_point(size = 3, color = "purple4", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  theme_classic(base_size = 15, base_family = "Times") +
  labs(
    x = "Castle ancestry proportion",
    y = expression("Inbreeding coefficient (" * italic(F) * ")")
  )

p_Fadmix

pdf("newadmixture_finbreeding_correlation.pdf",
    width = 6, height = 5,
    family = "Times",
    useDingbats = FALSE)

print(p_Fadmix)
dev.off()





#actual final graph for inbreeding 

library(dplyr)
library(ggplot2)

admixed_ids <- c("A03","B01","A02","A01","G01","E05","G03",
                 "D16","D27","D30","D19","D10","D12","D28")

fdat <- read.table("robustpopulation.het", header = TRUE) %>%
  mutate(
    Ho = (N_SITES - O.HOM.) / N_SITES,
    Group = case_when(
      INDV %in% admixed_ids ~ "Admixed",
      grepl("^C", INDV) ~ "Castle",
      grepl("^V", INDV) ~ "Wild",
      TRUE ~ "Northland"
    ),
    Group = factor(Group, levels = c("Wild", "Northland", "Castle", "Admixed"))
  )

# Summary table (inbreeding coefficient F)
f_summary <- fdat %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_F = mean(F, na.rm = TRUE),
    sd_F = sd(F, na.rm = TRUE),
    median_F = median(F, na.rm = TRUE),
    min_F = min(F, na.rm = TRUE),
    max_F = max(F, na.rm = TRUE),
    .groups = "drop"
  )

print(f_summary)

p_F <- ggplot(fdat, aes(x = Group, y = F, fill = Group)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.70,
    linewidth = 0.8
  ) +
  geom_jitter(
    aes(color = Group),
    width = 0.12,
    height = 0.01,
    size = 1.8,
    alpha = 0.75
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3.2,
    fill = "white",
    stroke = 0.8
  ) +
  scale_fill_manual(values = c(
    "Northland" = "#1f78b4",
    "Wild" = "grey",
    "Castle" = "#e31a1c",
    "Admixed" = "purple"
  )) +
  scale_color_manual(values = c(
    "Northland" = "#1f78b4",
    "Wild" = "grey30",
    "Castle" = "#e31a1c",
    "Admixed" = "purple4"
  )) +
  theme_classic(base_size = 15) +
  labs(
    x = NULL,
    y = expression("Inbreeding coefficient (" * italic(F) * ")")
  ) +
  theme(
    legend.position = "none",
    text = element_text(family = "Times New Roman"),
    axis.line = element_line(linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.9)
  )

print(p_F)

pdf("Inbreeding_boxplot_updated.pdf", width = 8, height = 6, family = "Times", useDingbats = FALSE)
print(p_F)
dev.off()




#heterozygosity correlation to castle ancestry 


f_admix <- fdat %>%
  left_join(admix, by = c("INDV" = "IID"))

admixed_Ho <- f_admix %>%
  filter(Castle > 0.02 & Castle < 0.95)

cor.test(admixed_Ho$Ho, admixed_Ho$Castle)
cor.test(admixed_Ho$Ho, admixed_Ho$Castle, method = "pearson")

p_Hoadmix <- ggplot(admixed_Ho, aes(x = Castle, y = Ho)) +
  geom_point(size = 3, color = "purple4", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  theme_classic(base_size = 15, base_family = "Times") +
  labs(
    x = "Castle ancestry proportion",
    y = expression("Observed heterozygosity (" * H[o] * ")")
  )
p_Hoadmix

pdf("d28admixture_heterozygosity_correlation.pdf",
    width = 6, height = 5,
    family = "Times",
    useDingbats = FALSE)

print(p_Hoadmix)
dev.off()


p_Hoadmix <- ggplot(admixed_Ho, aes(x = Castle, y = Ho)) +
  geom_jitter(
    width = 0.02,
    height = 0.005,
    size = 3,
    color = "purple4",
    alpha = 0.8
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  theme_classic(base_size = 15, base_family = "Times") +
  labs(
    x = "Castle ancestry proportion",
    y = expression("Observed heterozygosity (" * H[o] * ")")
  )

nrow(admixed_Ho)



model_Ho <- lm(Ho ~ Castle, data = admixed_Ho)
summary(model_Ho)
summary(model_Ho)$r.squared
summary(model_Ho)$coefficients[2,4]



# Inbreeding correlation to Castle ancestry

f_admix <- fdat %>%
  left_join(admix, by = c("INDV" = "IID"))

admixed_F <- f_admix %>%
  filter(Castle > 0.02 & Castle < 0.95)

cor.test(admixed_F$F, admixed_F$Castle, method = "pearson")

p_Fadmix <- ggplot(admixed_F, aes(x = Castle, y = F)) +
  geom_jitter(
    width = 0.02,
    height = 0.005,
    size = 3,
    color = "purple4",
    alpha = 0.8
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  theme_classic(base_size = 15, base_family = "Times") +
  labs(
    x = "Castle ancestry proportion",
    y = expression("Inbreeding coefficient (" * F * ")")
  )

p_Fadmix

pdf("d28admixture_inbreeding_correlation.pdf",
    width = 6, height = 5,
    family = "Times",
    useDingbats = FALSE)

print(p_Fadmix)
dev.off()


model_F <- lm(F ~ Castle, data = admixed_F)

summary(model_F)

# R-squared
summary(model_F)$r.squared

# slope
coef(model_F)[2]

# p-value
summary(model_F)$coefficients[2,4]

R² = summary(model_F)$r.squared
slope = coef(model_F)[2]
p = summary(model_F)$coefficients[2,4]
