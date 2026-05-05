library(ggplot2)
library(reshape2)

#fst_castle_legacy_noHighMiss.log:Weir and Cockerham weighted Fst estimate: 0.20923
#fst_castle_northland_noHighMiss.log:Weir and Cockerham weighted Fst estimate: 0.68722
#fst_northland_legacy_noHighMiss.log:Weir and Cockerham weighted Fst estimate: 0.52658


#Coromandel legacy
#fst_castle_legacy_noHighMiss.log:Weir and Cockerham weighted Fst estimate: 0.28152
#fst_castle_northland_noHighMiss.log:Weir and Cockerham weighted Fst estimate: 0.68722
#fst_northland_legacy_noHighMiss.log:Weir and Cockerham weighted Fst estimate: 0.60516

#Coro legacy with 80% missing removed (N=58)
#fst_castle_legacy.log:Weir and Cockerham weighted Fst estimate: 0.2717
#fst_castle_northland.log:Weir and Cockerham weighted Fst estimate: 0.67883
#fst_northland_legacy.log:Weir and Cockerham weighted Fst estimate: 0.60036

#no admixed individuals or data over 80% 12460 (5+5+39= total n49)
#fst_castle_legacyv2.log:Weir and Cockerham weighted Fst estimate: 0.2717
#fst_castle_northlandv2.log:Weir and Cockerham weighted Fst estimate: 0.70013
#fst_northland_legacyv2.log:Weir and Cockerham weighted Fst estimate: 0.61713

#after inclusion of bo1 and A04 (round them so 0.272, 0.7 and 0.613) 5+5+41= 
#fst_castle_legacy_updated.log:Weir and Cockerham weighted Fst estimate: 0.27196
#fst_castle_northland_updated.log:Weir and Cockerham weighted Fst estimate: 0.69554
#fst_northland_legacy_updated.log:Weir and Cockerham weighted Fst estimate: 0.61252


# Create symmetric matrix
fst_matrix <- matrix(
  c(
    0,        0.70013, 0.2717,
    0.70013,   0,      0.61313,
    0.2717,  0.61313, 0
  ),
  nrow = 3,
  byrow = TRUE
)

rownames(fst_matrix) <- c("Castle", "Northland", "Wild Coromandel")
colnames(fst_matrix) <- c("Castle", "Northland", "Wild Coromandel")

fst_long <- melt(fst_matrix)
colnames(fst_long) <- c("Pop1", "Pop2", "FST")

ggplot(fst_long, aes(Pop1, Pop2, fill = FST)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(FST == 0, "", round(FST, 3))),
            color = "white", size = 5) +
  scale_fill_gradient(low = "#deebf7", high = "#08306b") +
  labs(fill = expression(F[ST])) +   # <- HERE
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    text = element_text(family = "Times New Roman")
  ) +
  coord_fixed()

#trying to save as PDF
p <- ggplot(fst_long, aes(Pop1, Pop2, fill = FST)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(FST == 0, "", round(FST, 3))),
            color = "white", size = 5) +
  scale_fill_gradient(low = "#deebf7", high = "#08306b") +
  labs(fill = expression(F[ST])) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    text = element_text(family = "Times")   # <- use Times here
  ) +
  coord_fixed()

pdf("fst_heatmapv3.pdf", width = 6, height = 5, family = "Times", useDingbats = FALSE)
print(p)
dev.off()


fst_long$upper <- as.numeric(fst_long$Pop1) < as.numeric(fst_long$Pop2)

ggplot(subset(fst_long, upper), aes(Pop1, Pop2, fill = FST)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(FST, 3)), color = "white", size = 5) +
  scale_fill_gradient(low = "#deebf7", high = "#08306b") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    text = element_text(family = "Times New Roman")
  ) +
  coord_fixed()


ggsave("fst_plot.pdf", plot = p, device = cairo_pdf)







library(ggplot2)
library(reshape2)

fst_matrix <- matrix(
  c(
    0,        0.70013, 0.2717,
    0.70013,   0,      0.61313,
    0.2717,  0.61313, 0
  ),
  nrow = 3,
  byrow = TRUE
)

rownames(fst_matrix) <- c("Castle", "Northland", "Wild Coromandel")
colnames(fst_matrix) <- c("Castle", "Northland", "Wild Coromandel")

fst_long <- melt(fst_matrix)
colnames(fst_long) <- c("Pop1", "Pop2", "FST")

p <- ggplot(fst_long, aes(Pop1, Pop2, fill = FST)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(FST == 0, "", round(FST, 3))),
            color = "white", size = 5, family = "Times") +  # FIX HERE
  scale_fill_gradient(low = "#deebf7", high = "#08306b") +
  labs(fill = expression(F[ST])) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    text = element_text(family = "Times")
  ) +
  coord_fixed()

pdf("fst_heatmapv3.pdf", width = 6, height = 5, family = "Times", useDingbats = FALSE)
print(p)
dev.off()
