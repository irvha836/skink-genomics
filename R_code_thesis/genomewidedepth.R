library(dplyr)
library(ggplot2)

# ---- Load coverage windows ----
cov <- read.table(
  "coverage_100kb_ge5Mb.bed",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(cov) <- c("chrom", "start", "end", "depth")

# ---- Get contig lengths ----
lengths <- read.table(
  "contigs_ge5Mb.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

#thesis 
# Genome-wide depth plot from mosdepth windows (polished thesis version)
library(dplyr)
library(ggplot2)

# ---- Inputs ----
mosdepth_regions_gz <- "whitakers_q20.regions.bed.gz"   # mosdepth --by 100000 output
fai_path <- "whitakergenome.fasta.fai"                 # FASTA index for contig lengths

# OPTIONAL: exclude very short contigs to avoid window-edge artefacts at the tail end.
# 0 keeps everything; 1e5 = 100 kb; 5e5 = 500 kb (often a good compromise)
min_contig_len <- 5e5

# ---- Load mosdepth windows (chrom start end mean_depth) ----
cov <- read.table(gzfile(mosdepth_regions_gz), header = FALSE, stringsAsFactors = FALSE)
colnames(cov) <- c("chrom", "start", "end", "depth")

# ---- Load contig lengths from FASTA index (.fai) ----
fai <- read.table(fai_path, stringsAsFactors = FALSE)
colnames(fai)[1:2] <- c("chrom", "contig_length")

# ---- Pseudo-genome layout: order contigs largest -> smallest ----
contig_layout <- fai %>%
  arrange(desc(contig_length)) %>%
  mutate(genome_start = lag(cumsum(contig_length), default = 0)) %>%
  select(chrom, contig_length, genome_start)

# ---- Join depth windows to pseudo-genome coordinates ----
cov_genome <- cov %>%
  inner_join(contig_layout, by = "chrom") %>%
  mutate(
    genome_pos_bp = genome_start + start,
    genome_pos_gb = genome_pos_bp / 1e9
  ) %>%
  filter(is.finite(depth), depth >= 0)

# ---- Optional: drop ultra-short contigs (recommended for cleaner main figure) ----
if (!is.null(min_contig_len) && min_contig_len > 0) {
  cov_genome <- cov_genome %>% filter(contig_length >= min_contig_len)
}

# ---- X range (pseudo-genome length of what you are plotting) ----
total_gb <- sum(unique(cov_genome$contig_length)) / 1e9

# ---- Y cap for readability (points still exist above, just not shown) ----
y_cap <- quantile(cov_genome$depth, 0.99, na.rm = TRUE)

# ---- Use trimmed data for smoothing only (prevents spikes dragging LOESS) ----
cov_for_smooth <- cov_genome %>%
  filter(depth <= y_cap)

# ---- Global depth reference line ----
depth_global <- median(cov_genome$depth, na.rm = TRUE)

# ---- Plot ----
p <- ggplot(cov_genome, aes(x = genome_pos_gb, y = depth)) +
  geom_point(alpha = 0.08, size = 0.25) +
  geom_smooth(
    data = cov_for_smooth,
    method = "loess",
    span = 0.08,
    se = FALSE,
    linewidth = 1
  ) +
  geom_hline(yintercept = depth_global, linetype = "dashed") +
  coord_cartesian(xlim = c(0, total_gb),
                  ylim = c(15, y_cap)) +
  theme_classic(base_size = 15) +
  labs(
    x = "Genome position (Gb)",
    y = "Windowed depth (100 kb windows)"
  )

print(p)

# ---- Summary stats for Results ----
depth_mean   <- mean(cov_genome$depth, na.rm = TRUE)
depth_median <- median(cov_genome$depth, na.rm = TRUE)
depth_sd <- sd(cov_genome$depth, na.rm = TRUE)
depth_iqr <- IQR(cov_genome$depth, na.rm = TRUE)

# breadth of coverage in windows (0.5x–2x of median)
lo <- 0.5 * depth_median
hi <- 2.0 * depth_median
pct_in_range <- mean(cov_genome$depth >= lo & cov_genome$depth <= hi, na.rm = TRUE) * 100

depth_sd
depth_mean
depth_median
pct_in_range
depth_iqr

pct_zero <- mean(cov_genome$depth == 0, na.rm = TRUE) * 100
pct_zero

quantile(cov_genome$depth, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)





#investigating high coverage contigs 
cov_genome %>%
  arrange(desc(depth)) %>%
  select(chrom, start, end, depth) %>%
  head(20)



contig_depth %>%
  filter(chrom %in% c("contig_265_pilon_pilon",
                      "contig_2884_pilon_pilon",
                      "contig_2515_pilon_pilon")) 
fai %>%
  filter(chrom %in% c("contig_265_pilon_pilon",
                      "contig_2884_pilon_pilon",
                      "contig_2515_pilon_pilon"))


cov_genome %>%
  filter(chrom == "contig_265_pilon_pilon") %>%
  arrange(desc(depth)) %>%
  head(10)

cov_genome %>%
  filter(chrom == "contig_265_pilon_pilon") %>%
  summarise(mean_depth = mean(depth),
            median_depth = median(depth),
            max_depth = max(depth))






dist <- read.table("whitakers_q20.mosdepth.global.dist.txt",
                   header = FALSE,
                   stringsAsFactors = FALSE)

colnames(dist) <- c("depth", "bases", "prop")


str(dist)
head(dist)

library(ggplot2)

ggplot(dist, aes(x = depth, y = prop)) +
  geom_line(linewidth = 1) +
  theme_classic(base_size = 15) +
  xlab("Sequencing depth") +
  ylab("Proportion of genome")
#genome distribution depth 
dist <- read.table("whitakers_q20.mosdepth.global.dist.txt", header=FALSE)
colnames(dist) <- c("depth", "bases")

dist <- dist %>%
  mutate(prop = bases / sum(bases))

ggplot(dist, aes(x = depth, y = prop)) +
  geom_line(linewidth = 1) +
  theme_classic(base_size = 15) +
  xlab("Sequencing depth") +
  ylab("Proportion of genome")













library(dplyr)
library(ggplot2)

# 100 kb window depths
cov <- read.table(gzfile("whitakers_q20.regions.bed.gz"),
                  header = FALSE, stringsAsFactors = FALSE)
colnames(cov) <- c("chrom","start","end","depth")

# OPTIONAL: apply the same contig length filter as your main plot (recommended for consistency)
fai <- read.table("whitakergenome.fasta.fai", stringsAsFactors = FALSE)
colnames(fai)[1:2] <- c("chrom","contig_length")

min_contig_len <- 5e5  # set to 0 if you didn’t filter in your main plot
cov <- cov %>%
  inner_join(fai %>% select(chrom, contig_length), by="chrom") %>%
  filter(contig_length >= min_contig_len)

# histogram as % of genome (each window contributes equally; OK because all windows are 100 kb)
ggplot(cov, aes(x = depth)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                 bins = 80) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_classic(base_size = 15) +
  labs(x = "Mean depth (100 kb windows)",
       y = "Percentage of genome (%)")






ggplot(cov, aes(x = depth)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                 bins = 60) +
  coord_cartesian(xlim = c(15, 35)) +
  theme_classic(base_size = 15) +
  labs(x = "Mean depth (100 kb windows)",
       y = "Percentage of genome (%)")
















#plot v2 trial
library(dplyr)
library(ggplot2)

# ---- Inputs ----
mosdepth_regions_gz <- "whitakers_q20.regions.bed.gz"
fai_path <- "whitakergenome.fasta.fai"
min_contig_len <- 5e5   # keep contigs >= 500 kb

# ---- Load mosdepth windows ----
cov <- read.table(gzfile(mosdepth_regions_gz), header = FALSE, stringsAsFactors = FALSE)
colnames(cov) <- c("chrom", "start", "end", "depth")

# ---- Load contig lengths ----
fai <- read.table(fai_path, stringsAsFactors = FALSE)
colnames(fai)[1:2] <- c("chrom", "contig_length")

# ---- Pseudo-genome layout ----
contig_layout <- fai %>%
  arrange(desc(contig_length)) %>%
  mutate(genome_start = lag(cumsum(contig_length), default = 0)) %>%
  select(chrom, contig_length, genome_start)

# ---- Join windows to pseudo-genome coordinates ----
cov_genome <- cov %>%
  inner_join(contig_layout, by = "chrom") %>%
  mutate(
    genome_pos_bp = genome_start + start,
    genome_pos_gb = genome_pos_bp / 1e9
  ) %>%
  filter(is.finite(depth), depth >= 0)

# ---- Optional: remove very short contigs ----
if (!is.null(min_contig_len) && min_contig_len > 0) {
  cov_genome <- cov_genome %>% filter(contig_length >= min_contig_len)
}

# ---- Global depth line ----
depth_global <- median(cov_genome$depth, na.rm = TRUE)

# ---- Trimmed data only for smoothing ----
# keeps extreme spikes from pulling the smoother around
smooth_cap <- quantile(cov_genome$depth, 0.99, na.rm = TRUE)
cov_for_smooth <- cov_genome %>%
  filter(depth <= smooth_cap)

# =========================================================
# Main thesis figure: fixed axes to match text
# =========================================================
p_main <- ggplot(cov_genome, aes(x = genome_pos_gb, y = depth)) +
  geom_point(alpha = 0.08, size = 0.25) +
  geom_smooth(
    data = cov_for_smooth,
    method = "loess",
    span = 0.08,
    se = FALSE,
    linewidth = 1
  ) +
  geom_hline(yintercept = depth_global, linetype = "dashed") +
  coord_cartesian(
    xlim = c(0, 1.45),
    ylim = c(0, 50)
  ) +
  scale_x_continuous(
    breaks = seq(0, 1.4, by = 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 50, by = 10),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 15) +
  labs(
    x = "Genome position (Gb)",
    y = "Windowed depth (100 kb windows)"
  )

print(p_main)

ggsave(
  "genome_depth_main.pdf",
  plot = p_main,
  width = 9,
  height = 5,
  device = cairo_pdf
)



# Split datasets
cov_long <- cov_genome %>% filter(contig_length >= min_contig_len)
cov_short <- cov_genome %>% filter(contig_length < min_contig_len)

# Function to summarise
summarise_depth <- function(df) {
  data.frame(
    mean = mean(df$depth, na.rm = TRUE),
    median = median(df$depth, na.rm = TRUE),
    sd = sd(df$depth, na.rm = TRUE),
    cv = sd(df$depth, na.rm = TRUE) / mean(df$depth, na.rm = TRUE)
  )
}

long_stats <- summarise_depth(cov_long)
short_stats <- summarise_depth(cov_short)

long_stats
short_stats


cov_filtered <- cov_genome

if (!is.null(min_contig_len) && min_contig_len > 0) {
  cov_filtered <- cov_genome %>% 
    filter(contig_length >= min_contig_len)
}



depth_mean   <- mean(cov_filtered$depth, na.rm = TRUE)
depth_median <- median(cov_filtered$depth, na.rm = TRUE)
depth_ci     <- quantile(cov_filtered$depth, c(0.025, 0.975), na.rm = TRUE)

depth_mean
depth_median
depth_ci

x_max <- max(cov_filtered$genome_pos_gb, na.rm = TRUE)

p_filtered <- ggplot(cov_filtered, aes(x = genome_pos_gb, y = depth)) +
  geom_point(alpha = 0.08, size = 0.25) +
  geom_smooth(
    data = cov_filtered %>% 
      filter(depth <= quantile(depth, 0.99, na.rm = TRUE)),
    method = "loess",
    span = 0.08,
    se = FALSE,
    linewidth = 1
  ) +
  geom_hline(yintercept = depth_median, linetype = "dashed") +
  scale_x_continuous(
    limits = c(0, x_max),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 50),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 15) +
  labs(
    x = "Genome position (Gb)",
    y = "Windowed depth (100 kb windows)"
  )

print(p_filtered)



# =========================================================
# Appendix figure: full y-range so spikes are shown
# =========================================================
y_max_full <- ceiling(max(cov_genome$depth, na.rm = TRUE))

p_appendix <- ggplot(cov_genome, aes(x = genome_pos_gb, y = depth)) +
  geom_point(alpha = 0.08, size = 0.25) +
  geom_smooth(
    data = cov_for_smooth,
    method = "loess",
    span = 0.08,
    se = FALSE,
    linewidth = 1
  ) +
  geom_hline(yintercept = depth_global, linetype = "dashed") +
  coord_cartesian(
    xlim = c(0, 1.45),
    ylim = c(0, y_max_full)
  ) +
  scale_x_continuous(
    breaks = seq(0, 1.4, by = 0.2),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 15) +
  labs(
    x = "Genome position (Gb)",
    y = "Windowed depth (100 kb windows)"
  )

print(p_appendix)

ggsave(
  "genome_depth_appendix_fullrange.pdf",
  plot = p_appendix,
  width = 9,
  height = 5,
  device = cairo_pdf
)