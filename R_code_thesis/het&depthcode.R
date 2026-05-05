## Investigating heterozygosity and depth for Robust GBS data.

#import data from nesi 
#go session select working directory e.g downloads and load it in 
setwd ("~/Downloads")

# load in like this 
het <- read.table("outmindp3.het", header = TRUE)
depth <- read.table("outmindp3.idepth", header = TRUE)

# merge files by individuals as individuals are the same to produce
# merged data file
merged <- merge(het, depth, by = "INDV")


head(merged)
colnames(merged)

merged$H_obs <- 1 - (merged$O.HOM. / merged$N_SITES.x)

merged_lt50 <- subset(merged, MEAN_DEPTH < 50)

library(ggplot2)

ggplot(merged_lt50, aes(x = MEAN_DEPTH, y = H_obs)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Mean sequencing depth per individual (< 50)",
    y = "Observed heterozygosity",
    title = "Observed heterozygosity vs sequencing depth (majority of samples)"
  ) +
  theme_classic()

# now do same for mindp4 
het2 <- read.table("outmindp4.het", header = TRUE)
depth2 <- read.table("outmindp4.idepth", header = TRUE)

merged2 <- merge(het2, depth2, by = "INDV")

merged2$H_obs <- 1 - (merged2$O.HOM. / merged2$N_SITES.x)

merged2_lt50 <- subset(merged2, MEAN_DEPTH < 50)

ggplot(merged2_lt50, aes(x = MEAN_DEPTH, y = H_obs)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Mean sequencing depth per individual (< 50)",
    y = "Observed heterozygosity",
    title = "Observed heterozygosity vs sequencing depth mindp4 (majority of samples)"
  ) +
  theme_classic()


het60 <- read.table("out (2).het", header = TRUE)
depth60 <- read.table("out (2).idepth", header = TRUE)

merged3 <- merge(het60, depth60, by = "INDV")


merged3$H_obs <- 1 - (merged3$O.HOM. / merged3$N_SITES.x)

merged3_lt50 <- subset(merged3, MEAN_DEPTH < 50)

ggplot(merged3_lt50, aes(x = MEAN_DEPTH, y = H_obs)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Mean sequencing depth per individual (< 50)",
    y = "Observed heterozygosity",
    title = "Observed heterozygosity vs sequencing depth mindp3 (majority of samples)"
  ) +
  theme_classic()




# Vector of individuals to exclude
exclude_ids <- c("VUW01", "VUW02", "VUW03", "VUW04", "VUW05", "VUW06", "D27", "D19", "C04", "D06", "D15", "D17", "D28", "E02", "F11")

# Subset the data frame to exclude those individuals
filtered <- subset(merged, !(INDV %in% exclude_ids))

# Plot heterozygosity (F) vs sequencing depth
plot(filtered$MEAN_DEPTH, filtered$F,
     xlab = "Mean Depth", ylab = "Inbreeding Coefficient (F)",
     main = "Heterozygosity vs Depth (High-depth individuals removed)")

#remove individuals with high depth to tease apart relationship on small scale
# (VUW01, VUW02, VUW03, VUW04, VUW05, VUW06, D27, D19)
#plot merged file mean depth against heterozygosity (F)
plot(merged$MEAN_DEPTH, merged$F, 
     xlab = "Mean Depth", ylab = "Inbreeding Coefficient (F)",
     main = "Heterozygosity vs Depth")




#or can run a ggplot but first have to laod ggplot
library(ggplot2)

ggplot(filtered, aes(x = MEAN_DEPTH, y = F)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Mean Depth", y = "Inbreeding Coefficient (F)",
       title = "Heterozygosity vs Sequencing Depth (Excluding High-Depth Individuals)") +
  theme_minimal()

merged4 <- merge(hetmindp4, depthmindp4, by = "INDV")
merged5 <- merge(hetmindp5, depthmindp5, by = "INDV")
merged6 <- merge(hetmindp6, depthmindp6, by = "INDV")

filtered4 <- subset(merged4, !(INDV %in% exclude_ids))
filtered5 <- subset(merged5, !(INDV %in% exclude_ids))
filtered6 <- subset(merged6, !(INDV %in% exclude_ids))


ggplot(filtered4, aes(x = MEAN_DEPTH, y = F)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Mean Depth", y = "Inbreeding Coefficient (F)",
       title = "Heterozygosity vs Sequencing Depth (Excluding High-Depth Individuals)") +
  theme_minimal()


ggplot(filtered5, aes(x = MEAN_DEPTH, y = F)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Mean Depth", y = "Inbreeding Coefficient (F)",
       title = "Heterozygosity vs Sequencing Depth (Excluding High-Depth Individuals)") +
  theme_minimal()


ggplot(filtered6, aes(x = MEAN_DEPTH, y = F)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Mean Depth", y = "Inbreeding Coefficient (F)",
       title = "Heterozygosity vs Sequencing Depth (Excluding High-Depth Individuals)") +
  theme_minimal()
# can see as men depth increases homozygosity decreases therefore individuals
# with higher sequencing depth display more heterozygosity and sequencing depth
# is not independent from heterozygosity so low sequencing depth is likely not ideal 
# as it underestimates heterozygosity - maybe i need to be more ruthless in my Filtering