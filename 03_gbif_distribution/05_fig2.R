# Step 05: Fig. 2 -- (a) global map of thinned occurrences coloured by anthocyanin
# status, (b) absolute latitude by anthocyanin status.
#
# Input : data/SupplementaryData02_GBIF_thinned_occurrence_records.csv
# Output: output/Fig2.pdf
# Run from this directory:  Rscript 05_fig2.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(maps)
  library(ggbeeswarm)
  library(patchwork)
})

records <- readr::read_csv("data/SupplementaryData02_GBIF_thinned_occurrence_records.csv",
                           show_col_types = FALSE)

world_map <- ggplot2::map_data("world")

antho_colors <- c("Antho" = "#9076b4", "NonAntho" = "#719d45")
antho_labels <- c("Antho" = "Anthocyanin accumulating",
                  "NonAntho" = "Anthocyanin non-accumulating")

records <- records |>
  mutate(antho = factor(antho, levels = c("Antho", "NonAntho")),
         abs_latitude = abs(decimalLatitude)) |>
  arrange(antho)

theme_fig2 <- theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, color = "black"),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(4, "mm"),
    plot.tag = element_text(size = 11, face = "bold"),
    plot.tag.position = c(0.01, 0.98),
    plot.margin = margin(3, 3, 3, 3)
  )

p_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "grey88", color = "white", linewidth = 0.08) +
  geom_point(data = records,
             aes(x = decimalLongitude, y = decimalLatitude, color = antho),
             size = 0.22, alpha = 0.55, stroke = 0) +
  scale_color_manual(values = antho_colors, labels = antho_labels, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_fig2 +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.25),
    panel.grid.major = element_line(color = "grey82", linewidth = 0.18),
    panel.grid.minor = element_blank()
  )

p_box <- ggplot(records, aes(x = antho, y = abs_latitude, fill = antho, color = antho)) +
  geom_boxplot(width = 0.52, alpha = 0.45, outlier.shape = NA, linewidth = 0.35) +
  ggbeeswarm::geom_quasirandom(size = 0.22, alpha = 0.22, width = 0.18, stroke = 0) +
  scale_fill_manual(values = antho_colors, labels = antho_labels, drop = FALSE) +
  scale_color_manual(values = antho_colors, labels = antho_labels, drop = FALSE) +
  scale_x_discrete(labels = c("Antho" = "Accum.", "NonAntho" = "Non-accum.")) +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 30),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Absolute latitude") +
  theme_fig2 +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.25)
  )

fig2 <- (p_map + p_box + plot_layout(widths = c(2.9, 1))) +
  plot_annotation(tag_levels = "A") &
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("output/Fig2.pdf", fig2, width = 7.2, height = 3.45, units = "in", device = "pdf")
