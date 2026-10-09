# Step 03: Fig. 4c - light-response curves of Y(II), ETR and NPQ (mean of the
#          three plant samples per species x growth light at each actinic PPFD).
#
# Input : output/PAM_light_curves_processed.csv     (01)
# Output: output/Fig4c.pdf
#
# Plotting code and styles taken from the submission script that assembled Fig. 4;
# the final panel arrangement and species abbreviations in the published figure
# were set in the layout software.
# Author: Natsu Katayama
# Run from this directory:  Rscript 03_fig4c_light_response_curves.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

OUT <- "output"

theme_base <- theme_bw(base_size = 7.6) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.6),
    plot.margin = margin(2, 4, 2, 4)
  )

species_labels   <- c("Sp7498" = "Spirodela polyrhiza", "Wa8730" = "Wolffia australiana")
species_colors   <- c("Sp7498" = "#8c77af", "Wa8730" = "#73a244")
growth_labels    <- c("PFD_50" = "LL", "PFD_1000" = "HL")
growth_shapes    <- c("PFD_50" = 16, "PFD_1000" = 17)
growth_linetypes <- c("PFD_50" = "solid", "PFD_1000" = "dashed")

lc_summary <- read.csv(file.path(OUT, "PAM_light_curves_processed.csv")) |>
  mutate(
    SpID = factor(SpID, levels = names(species_labels)),
    pfd = factor(pfd, levels = names(growth_labels))
  ) |>
  group_by(SpID, pfd, PPFD) |>
  summarise(
    `Y(II)` = mean(Y.II., na.rm = TRUE),
    ETR = mean(ETR2, na.rm = TRUE),
    NPQ = mean(NPQ, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(`Y(II)`, ETR, NPQ), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = c("Y(II)", "ETR", "NPQ")))

# The last light step (PPFD 3005.5) lies outside the plotted range (0-2000).
pC <- ggplot(lc_summary, aes(x = PPFD, y = value, group = interaction(SpID, pfd))) +
  geom_line(aes(color = SpID, linetype = pfd), linewidth = 0.62, alpha = 0.95) +
  geom_point(aes(color = SpID, shape = pfd), size = 1.70, stroke = 0, alpha = 0.95) +
  facet_wrap(~ parameter, nrow = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0, 2000)) +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  scale_shape_manual(values = growth_shapes, labels = growth_labels, name = "Growth light") +
  scale_linetype_manual(values = growth_linetypes, labels = growth_labels, name = "Growth light") +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 16, linetype = 0, size = 2.3)),
    shape = guide_legend(order = 2, override.aes = list(color = "black", linetype = 0, size = 2.1)),
    linetype = "none"
  ) +
  labs(x = expression("Actinic PPFD (" * mu * "mol photons " * m^-2 * " " * s^-1 * ")"), y = NULL) +
  theme_base +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.3),
    strip.text = element_text(size = 7.4),
    legend.text = element_text(size = 6.6),
    axis.title = element_text(size = 8.0),
    axis.text = element_text(size = 7.1, colour = "black")
  )

ggsave(file.path(OUT, "Fig4c.pdf"), pC, width = 7.2, height = 3.0, units = "in", device = "pdf")
message("Wrote: ", file.path(OUT, "Fig4c.pdf"))
