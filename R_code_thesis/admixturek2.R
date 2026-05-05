library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Downloads")

fam_file <- "robust_captive_pruned.fam"
q2_file  <- "robust_captive_pruned.2.Q"

# ----------------------------
# Read .fam and define a priori population labels
# ----------------------------
fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fam) <- c("FID","IID","PAT","MAT","SEX","PHENO")

fam <- fam %>%
  mutate(Population = ifelse(grepl("^C", IID), "Castle", "Northland"))

# ----------------------------
# Read Q2 (wide)
# ----------------------------
q2 <- read.table(q2_file, header = FALSE, stringsAsFactors = FALSE)
colnames(q2) <- paste0("Cluster", seq_len(ncol(q2)))

q2 <- bind_cols(fam %>% select(IID, Population), q2)


# ----------------------------
# Relabel clusters so "Castle" and "Northland" components are consistent
# ----------------------------
relabel_clusters_K2 <- function(df){
  cl <- grep("^Cluster", names(df), value = TRUE)
  
  castle_means <- df %>%
    filter(Population == "Castle") %>%
    summarise(across(all_of(cl), ~mean(.x, na.rm = TRUE)))
  
  castle_cl <- cl[which.max(as.numeric(castle_means[1, ]))]
  north_cl  <- setdiff(cl, castle_cl)[1]
  
  out <- df
  names(out)[names(out) == castle_cl] <- "Castle"
  names(out)[names(out) == north_cl]  <- "Northland"
  out
}

q2r <- relabel_clusters_K2(q2)

q2r %>% select(IID, Population, Castle, Northland)
# ----------------------------
# Order individuals: Castle first, then Northland
# Within each group, sort by Castle ancestry (from K=2) for a clean gradient
# ----------------------------
order_ids <- q2r %>%
  mutate(Castle = ifelse(is.na(Castle), 0, Castle)) %>%
  arrange(Population, desc(Castle), IID) %>%
  pull(IID)

n_castle  <- sum(fam$Population == "Castle")
divider_x <- n_castle + 0.5

# ----------------------------
# Long format for plotting
# ----------------------------
admix_long <- q2r %>%
  select(IID, Population, Castle, Northland) %>%
  pivot_longer(
    cols = c(Castle, Northland),
    names_to = "Component",
    values_to = "Ancestry"
  ) %>%
  mutate(
    IID       = factor(IID, levels = order_ids),
    Component = factor(Component, levels = c("Castle","Northland"))
  )

# ----------------------------
# Colours (match your original)
# ----------------------------
cols <- c(Castle = "#d73027", Northland = "#4575b4")

# ----------------------------
# Plot (styled like your original K2/K3 figure)
# ----------------------------
p <- ggplot(admix_long, aes(x = IID, y = Ancestry, fill = Component)) +
  geom_col(width = 1, color = "white", linewidth = 0.1) +
  
  # Divider between groups
  annotate("segment", x = divider_x, xend = divider_x, y = 0, yend = 1,
           linewidth = 1.2, color = "black") +
  
  # Group labels above bars
  annotate("text", x = n_castle/2, y = 1.06, label = "Castle",
           fontface = "bold", size = 5) +
  annotate("text",
           x = n_castle + (length(order_ids)-n_castle)/2, y = 1.06,
           label = "Northland", fontface = "bold", size = 5) +
  
  scale_fill_manual(values = cols, drop = FALSE) +
  
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1.08),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.margin  = margin(10, 10, 25, 10)
  ) +
  labs(x = NULL, y = "Ancestry proportion", fill = "Ancestry")

print(p)

# --- make sure your plot uses the macOS Times mapping
p <- p + theme_classic(base_size = 14, base_family = "Times")

# --- save to PDF (no Cairo/XQuartz needed)
pdf("admixture_K2v2.pdf", width = 10, height = 4, family = "Times", useDingbats = FALSE)
print(p)
dev.off()

p <- p +
  theme_classic(base_size = 14, base_family = "Times") +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.margin  = margin(10, 10, 25, 10)
  )



ggsave("admixture_K2_pub.png", p, width = 12, height = 3.8, dpi = 600)
ggsave("admixture_K2_pub.pdf", p, width = 12, height = 3.8)
# If you need cairo:
ggsave("admixture_K2_pub.pdf", p, width = 12, height = 3.8, device = cairo_pdf)


ggsave(
  "fst_heatmap.pdf",
  width = 6,
  height = 5,
  device = cairo_pdf
)

head Q
