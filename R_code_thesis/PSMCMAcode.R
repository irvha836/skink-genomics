#PSMC code for thesis
# used code from whooping crane paper as inspo
# available here https://github.com/claudefa/WhoopingCrane_genomics/blob/main/Rscripts/Figure3.R
library(ggplot2)
library(scales)
library(psych)  # for harmonic.mean

# ---- Functions to read and scale PSMC ----
read_psmc <- function(x){
  out <- scan(x, what = "", sep = "\n", quiet = TRUE)
  niters <- as.integer(gsub("^.*n_iterations:|,.*$", "", out[14]))
  parapattern <- gsub("^.*pattern:|,.*$", "", out[13])
  nintervs <- eval(parse(text = parapattern))
  decoding <- as.integer(gsub("^.*is_decoding:", "", out[15]))
  res <- list()
  res$niters <- niters
  res$n <- nintervs
  # extract RS lines
  s <- grep("^RS", out, value = TRUE)
  s <- gsub("^RS\t", "", s)
  s <- strsplit(s, "\t")
  s <- matrix(as.numeric(unlist(s)), ncol = 6, byrow = TRUE)
  colnames(s) <- c("k","t_k","lambda_k","pi_k","sum_A_kl","A_kk")
  iter <- rep(0:res$niters, each = res$n)
  s <- cbind(s, iter = iter)
  res$RS <- s
  # extract theta
  s <- grep("^TR", out, value = TRUE)
  s <- gsub("^TR\t", "", s)
  s <- as.numeric(unlist(strsplit(s, "\t")))
  res$theta0 <- s[1]
  res$nintervs <- nintervs
  res
}

getXYplot.psmc <- function(x, mutation.rate, g, bin.size){
  RS <- x$RS[x$RS[, "iter"] == x$niters, ]
  N0 <- x$theta0 / (4*mutation.rate*bin.size)
  xx <- 2 * N0 * RS[, "t_k"] * g   # years
  yy <- N0 * RS[, "lambda_k"]       # Ne
  list(xx = xx, yy = yy)
}

# ---- Load main PSMC ----
psmc_skink <- read_psmc("skink.psmc")

# ---- Get X/Y for plotting ----
mutation_rate <- 1.17e-8
generation_time <- 6.5
bin_size <- 100  # you can adjust
xy_main <- getXYplot.psmc(psmc_skink, mutation_rate, generation_time, bin_size)
df_main <- data.frame(time = xy_main$xx, Ne = xy_main$yy)

# ---- Load bootstrap PSMCs ----
bootstrap_files <- list.files("rounds2", pattern="\\.psmc$", full.names = TRUE)

getXYplot_boot <- function(psmc_file, mutation.rate, g, bin.size, replicate_id){
  psmc <- read_psmc(psmc_file)
  xy <- getXYplot.psmc(psmc, mutation.rate, g, bin.size)
  data.frame(time = xy$xx, Ne = xy$yy, replicate = replicate_id)
}

boot_list <- lapply(seq_along(bootstrap_files), function(i){
  getXYplot_boot(bootstrap_files[i], mutation_rate, generation_time, bin_size, i)
})
df_boot <- do.call(rbind, boot_list)

# =========================
# PLOT IN Ma (no decimals on axis labels)
# =========================
# Convert years -> Ma for plotting
df_main$time_Ma <- df_main$time / 1e6
df_boot$time_Ma <- df_boot$time / 1e6

# Axis breaks you asked for (1–30 Ma)
ma_breaks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30)

ggplot() +
  geom_step(data = df_boot, aes(x = time_Ma, y = Ne, group = replicate),
            color="grey70", alpha=0.4) +
  geom_step(data = df_main, aes(x = time_Ma, y = Ne),
            color="red", size=1.2) +
  scale_x_log10(
    breaks = ma_breaks,
    labels = ma_breaks
  ) +
  scale_y_continuous(labels = scales::comma) +
  xlab("Millions of Years Ago (Ma)") +
  ylab("Effective Population Size (Ne)") +
  theme_classic()

# PDF for publication
pdf("PSMC_Skink_plot_Ma.pdf", width = 6, height = 4)  # width & height in inches

ggplot() +
  geom_step(data = df_boot, aes(x = time_Ma, y = Ne, group = replicate),
            color="grey70", alpha=0.4) +
  geom_step(data = df_main, aes(x = time_Ma, y = Ne),
            color="red", size=0.8) +
  scale_x_log10(
    breaks = ma_breaks,
    labels = ma_breaks
  ) +
  scale_y_continuous(labels = scales::comma) +
  xlab("Millions of Years Ago (Ma)") +
  ylab("Effective Population Size (Ne)") +
  theme_classic(base_size = 12)

dev.off()

tiff("PSMC_Skink_plot_Ma.tiff", width=6, height=4, units="in", res=600)
ggplot() +
  geom_step(data = df_boot, aes(x = time_Ma, y = Ne, group = replicate),
            color="grey70", alpha=0.4) +
  geom_step(data = df_main, aes(x = time_Ma, y = Ne),
            color="red", size=0.8) +
  scale_x_log10(
    breaks = ma_breaks,
    labels = ma_breaks
  ) +
  scale_y_continuous(labels = scales::comma) +
  xlab("Millions of Years Ago (Ma)") +
  ylab("Effective Population Size (Ne)") +
  theme_classic(base_size = 12)
dev.off()

# =========================
# Summary stats
# =========================
harmonic_Ne <- harmonic.mean(df_main$Ne)
harmonic_Ne

harmonic_Ne_1M <- harmonic.mean(df_main$Ne[df_main$time <= 1e6])
harmonic_Ne_1M

max_Ne <- max(df_main$Ne)
max_Ne

max_point <- df_main[which.max(df_main$Ne), ]
max_point

min_Ne <- min(df_main$Ne)
min_Ne

min_Ne_recent <- min(df_main$Ne[df_main$time <= 1e6])
min_Ne_recent

range(df_main$time) / 1e6   # min/max time in millions of years (Ma)
max(df_main$time) / 1000    # max time in kya

harmonic.mean(df_main$Ne[df_main$time <= 10e6])
range(df_main$time[df_main$time <= 100000])

present_like <- df_main[which.min(df_main$time), ]
present_like

Ne_most_recent <- df_main$Ne[which.min(df_main$time)]
Ne_most_recent

Ne_most_recent <- df_main %>%
  arrange(time) %>%
  slice(1)
Ne_most_recent

min_point <- df_main[which.min(df_main$Ne), ]
min_point
