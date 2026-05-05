library(dplyr)
library(ggplot2)


## ----------------------------
## 0) REQUIRED INPUTS
## ----------------------------
## het must contain: chrom, window_start, calls, hets
## contig_layout must contain: chrom, genome_start  (both in bp)

het <- read.table(
  "genome_wide_het_windows.txt.gz",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Rename columns (adjust if necessary)
colnames(het) <- c("chrom", "window_start", "sites_total", "calls", "hets")

fai <- read.table("whitakergenome_hardmasked.fasta.fai", stringsAsFactors = FALSE)
colnames(fai) <- c("chrom", "contig_length", "offset", "linebases", "linewidth")

contig_layout <- fai %>%
  arrange(desc(contig_length)) %>%
  mutate(
    genome_start = lag(cumsum(contig_length), default = 0)
  ) %>%
  select(chrom, genome_start)


## ----------------------------
## 1) Window-level QC + map to genome coordinates
## ----------------------------
het_genome <- het %>%
  transmute(
    chrom = chrom,
    window_start = window_start,
    calls = as.numeric(calls),
    hets  = as.numeric(hets)
  ) %>%
  filter(
    is.finite(calls), is.finite(hets),
    calls > 0, hets >= 0,
    calls >= 1000               # window-level minimum callability
  ) %>%
  inner_join(contig_layout %>% select(chrom, genome_start), by = "chrom") %>%
  mutate(genome_pos = genome_start + window_start)  # bp, cumulative position

## genome-wide pi (weighted by callability)
genome_pi <- sum(het_genome$hets) / sum(het_genome$calls)

## ----------------------------
## 2) Bin into 10 Mb and compute pi per bin
## ----------------------------
bin_size <- 1e7  # 10 Mb

het_binned <- het_genome %>%
  mutate(bin = floor(genome_pos / bin_size) * bin_size) %>%
  group_by(bin) %>%
  summarise(
    total_hets  = sum(hets),
    total_calls = sum(calls),
    pi = total_hets / total_calls,
    .groups = "drop"
  ) %>%
  arrange(bin) %>%
  mutate(x_gb = bin / 1e9)

## ----------------------------
## 3) OPTIONAL: mask unreliable bins + cap rare spikes
##    (keeps the plot stable if you have low-call tail bins)
## ----------------------------
min_calls <- 0.2 * median(het_binned$total_calls, na.rm = TRUE)

pi_cap <- quantile(
  het_binned$pi[het_binned$total_calls >= min_calls],
  probs = 0.995,
  na.rm = TRUE
)

plot_df <- het_binned %>%
  mutate(
    pi_plot = ifelse(total_calls >= min_calls & pi <= pi_cap, pi, NA_real_)
  )

## ----------------------------
## 4) Plot (bars + dashed genome mean)
## ----------------------------
ggplot(plot_df, aes(x = x_gb, y = pi_plot)) +
  geom_col(width = 0.008, na.rm = TRUE) +
  geom_hline(yintercept = genome_pi, linetype = "dashed", linewidth = 0.8) +
  theme_classic(base_size = 15) +
  scale_x_continuous(
    breaks = seq(0, max(plot_df$x_gb, na.rm = TRUE), by = 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Genome position (Gb)",
    y = expression("Nucleotide diversity (" * pi * ", 10 Mb bins)")
  )



ggplot(het_binned_clean, aes(x = x_gb, y = pi_plot)) +
  geom_col(width = 0.008, fill = "#2C7FB8", na.rm = TRUE) +
  geom_hline(yintercept = genome_pi,
             linetype = "dashed", linewidth = 0.8, color = "grey35") +
  scale_x_continuous(
    breaks = seq(0, max(het_binned_clean$x_gb, na.rm = TRUE), by = 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    plot.margin = margin(10, 12, 10, 10)
  ) +
  labs(
    x = "Genome position (Gb)",
    y = expression("Nucleotide diversity (" * pi * ", 10 Mb bins)")
  )



bin_size <- 1e7

hb <- plot_df %>% arrange(bin)   # or het_binned, but use the one you PLOT

hb %>%
  mutate(gap_bp = bin - lag(bin)) %>%
  filter(!is.na(gap_bp) & gap_bp != bin_size) %>%
  select(bin, x_gb, gap_bp) %>%
  head(20)





plot_df %>%
  filter(x_gb > 0.38, x_gb < 0.50) %>%
  summarise(
    n_bins = n(),
    n_na = sum(is.na(pi_plot)),
    min_calls = min(total_calls),
    median_calls = median(total_calls)
  )

plot_df %>%
  filter(x_gb > 0.38, x_gb < 0.50, is.na(pi_plot)) %>%
  select(bin, x_gb, total_calls, pi)



het_binned %>%
  filter(x_gb > 0.38 & x_gb < 0.42)





plot_df <- het_binned %>%
  mutate(
    is_outlier = pi > pi_cap,
    is_lowcall = total_calls < min_calls,
    is_flagged = is_outlier | is_lowcall
  )

ggplot(plot_df, aes(x = x_gb, y = pi, fill = is_flagged)) +
  geom_col(width = 0.01) +
  geom_hline(yintercept = genome_pi, linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = c("FALSE" = "#2C7FB8", "TRUE" = "grey80")) +
  guides(fill = "none") +
  theme_classic(base_size = 15) +
  scale_x_continuous(
    breaks = seq(0, max(plot_df$x_gb, na.rm = TRUE), by = 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Genome position (Gb)",
    y = expression("Nucleotide diversity (" * pi * ", 10 Mb bins)")
  )


max_x <- 1.2  # choose cut before the weird tail begins
plot_main <- plot_df %>% filter(x_gb <= max_x)

ggplot(plot_main, aes(x = x_gb, y = pi, fill = is_flagged)) +
  geom_col(width = 0.01) +
  geom_hline(yintercept = genome_pi, linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = c("FALSE" = "#2C7FB8", "TRUE" = "grey80")) +
  guides(fill = "none") +
  theme_classic(base_size = 15) +
  labs(x = "Genome position (Gb)",
       y = expression("Nucleotide diversity (" * pi * ", 10 Mb bins)"))

min_contig_included <- het_genome %>%
  left_join(fai %>% select(chrom, contig_length), by = "chrom") %>%
  summarise(min_len = min(contig_length, na.rm = TRUE)) %>%
  pull(min_len)

min_contig_included



mean_bin_pi <- mean(het_binned$pi, na.rm = TRUE)
sd_bin_pi   <- sd(het_binned$pi, na.rm = TRUE)
median_bin_pi <- median(het_binned$pi, na.rm = TRUE)

range_bin_pi <- range(het_binned$pi, na.rm = TRUE)




range_bin_pi 
median_bin_pi
sd_bin_pi
mean_bin_pi
sd_masked <- sd(plot_df$pi_plot, na.rm = TRUE)



plot_df <- het_binned %>%
  mutate(
    is_outlier = pi > pi_cap,
    is_lowcall = total_calls < min_calls,
    is_flagged = is_outlier | is_lowcall
  )
sd_masked <- sd(plot_df$pi[!plot_df$is_flagged], na.rm = TRUE)
mean_masked <- mean(plot_df$pi[!plot_df$is_flagged], na.rm = TRUE)
median_masked <- median(plot_df$pi[!plot_df$is_flagged], na.rm = TRUE)
range_masked <- range(plot_df$pi[!plot_df$is_flagged], na.rm = TRUE)





sd_masked 
mean_masked 
median_masked 
range_masked



quantile(plot_df$pi[!plot_df$is_flagged], probs = c(0.025, 0.975), na.rm = TRUE)

genome_pi_masked <- sum(plot_df$total_hets[!plot_df$is_flagged]) /
  sum(plot_df$total_calls[!plot_df$is_flagged])
