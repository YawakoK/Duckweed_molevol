# Step 06: Supplementary Fig. 1 -- Lemna only: (page 1) map of thinned occurrences
# coloured by anthocyanin status, (page 2) absolute latitude by anthocyanin status.
#
# Input : data/SupplementaryData02_GBIF_thinned_occurrence_records.csv
# Output: output/SupplementaryFig1.pdf (2 pages, A4)
#         ../environment/R_sessionInfo_03_gbif_distribution.txt
# Run from this directory:  Rscript 06_suppfig1.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(maps)
  library(ggbeeswarm)
})

lemna <- read.csv("data/SupplementaryData02_GBIF_thinned_occurrence_records.csv") |>
  filter(genus_cor == "Lemna") |>
  mutate(antho = factor(antho, levels = c("Antho", "NonAntho")))

world_map <- map_data("world")

antho_colors <- c("Antho" = "#9076b4", "NonAntho" = "#719d45")
antho_labels <- c("Antho" = "Anthocyanin present", "NonAntho" = "Anthocyanin absent")

p_lemna_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white", linewidth = 0.1) +
  geom_point(data = arrange(lemna, antho),
             aes(x = decimalLongitude, y = decimalLatitude, color = antho),
             size = 0.4, alpha = 0.6) +
  scale_color_manual(values = antho_colors, labels = antho_labels) +
  theme_bw() +
  ggtitle("Lemna distribution by anthocyanin status") +
  labs(x = "Longitude", y = "Latitude", color = NULL) +
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 85)) +
  theme(panel.grid = element_blank())

p_lemna_box <- ggplot(lemna, aes(x = antho, y = abs(decimalLatitude), fill = antho)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_quasirandom(aes(color = antho), size = 0.4, alpha = 0.25, width = 0.2) +
  scale_fill_manual(values = antho_colors, labels = antho_labels) +
  scale_color_manual(values = antho_colors, labels = antho_labels) +
  scale_x_discrete(labels = antho_labels) +
  theme_bw() +
  ggtitle("Absolute latitude of Lemna by anthocyanin status") +
  labs(x = "Anthocyanin status", y = "Absolute latitude") +
  theme(panel.grid = element_blank(), legend.position = "none")

pdf("output/SupplementaryFig1.pdf", width = 8.3, height = 11.7)
print(p_lemna_map)
print(p_lemna_box)
invisible(dev.off())

writeLines(capture.output(sessionInfo()), "../environment/R_sessionInfo_03_gbif_distribution.txt")
