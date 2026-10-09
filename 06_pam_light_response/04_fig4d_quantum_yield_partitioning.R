# Step 04: Fig. 4d - Y(II), Y(NPQ) and Y(NO) at AL50 and AL1000 for LL- and HL-grown
#          plants (mean +/- SE, n = 3), with S. polyrhiza vs W. australiana brackets
#          (Sidak-adjusted contrasts from Supplementary Data 16).
#
# Input : output/SupplementaryData15_PAM_quantum_yield_summary.csv          (02)
#         output/SupplementaryData16_PAM_species_pairwise_contrasts_by_condition.csv (02)
# Output: output/Fig4d.pdf
#
# Plotting code and styles taken from the submission script that assembled Fig. 4;
# P-value typography in the published figure was edited in the layout software.
# Author: Natsu Katayama
# Run from this directory:  Rscript 04_fig4d_quantum_yield_partitioning.R
suppressPackageStartupMessages({
  library(dplyr)
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

ytype_levels <- c("Y(II)", "Y(NPQ)", "Y(NO)")

quant_summary <- read.csv(file.path(OUT, "SupplementaryData15_PAM_quantum_yield_summary.csv"),
                          check.names = FALSE) |>
  mutate(
    species = factor(species, levels = c("Sp", "Wa")),
    growth_light = factor(growth_light, levels = c("LL", "HL")),
    actinic_light = factor(actinic_light, levels = c("AL50", "AL1000")),
    Ytype = factor(Ytype, levels = ytype_levels),
    growth_x = if_else(growth_light == "LL", 1.00, 1.68)
  )

spwa_stats <- read.csv(file.path(OUT, "SupplementaryData16_PAM_species_pairwise_contrasts_by_condition.csv"),
                       check.names = FALSE) |>
  mutate(
    growth_light = factor(growth_light, levels = c("LL", "HL")),
    actinic_light = factor(actinic_light, levels = c("AL50", "AL1000")),
    Ytype = factor(Ytype, levels = ytype_levels),
    label = paste0(p_label, "\n", case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "n.s."
    )),
    x_bracket = if_else(growth_light == "LL", 0.86, 1.82)
  )

ypos <- quant_summary |>
  group_by(actinic_light, Ytype, growth_light) |>
  summarise(
    y_start = max(mean[species == "Sp"] + se[species == "Sp"], na.rm = TRUE),
    y_end = max(mean[species == "Wa"] + se[species == "Wa"], na.rm = TRUE),
    .groups = "drop"
  )

spwa_stats <- spwa_stats |>
  select(actinic_light, Ytype, growth_light, label, x_bracket) |>
  left_join(ypos, by = c("actinic_light", "Ytype", "growth_light")) |>
  mutate(
    y_start = pmax(y_start, 0.02),
    y_end = pmax(y_end, 0.02),
    y_mid = (y_start + y_end) / 2,
    text_x = x_bracket + if_else(growth_light == "LL", -0.035, 0.035),
    text_hjust = if_else(growth_light == "LL", 1, 0)
  )

pD <- ggplot(quant_summary, aes(x = growth_x, y = mean, color = species, shape = species, group = species)) +
  geom_line(linewidth = 0.70, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.16, linewidth = 0.62) +
  geom_point(size = 1.90) +
  geom_segment(
    data = spwa_stats,
    aes(x = x_bracket, xend = x_bracket, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "grey25"
  ) +
  geom_text(
    data = spwa_stats,
    aes(x = text_x, y = y_mid, label = label, hjust = text_hjust),
    inherit.aes = FALSE,
    size = 1.90,
    lineheight = 0.85,
    color = "grey20"
  ) +
  facet_grid(actinic_light ~ Ytype) +
  scale_x_continuous(
    breaks = c(1.00, 1.68),
    labels = c("LL", "HL"),
    limits = c(0.68, 2.02),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_manual(values = c("Sp" = "#8c77af", "Wa" = "#73a244"),
                     labels = c("Sp" = "Spirodela polyrhiza", "Wa" = "Wolffia australiana"),
                     name = "Species") +
  scale_shape_manual(values = c("Sp" = 15, "Wa" = 16),
                     labels = c("Sp" = "Spirodela polyrhiza", "Wa" = "Wolffia australiana"),
                     name = "Species") +
  coord_cartesian(ylim = c(0, 1.08), clip = "off") +
  labs(x = "Growth condition", y = "Quantum yield") +
  theme_base +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 7.4),
    legend.position = "top",
    legend.text = element_text(size = 6.8),
    axis.title = element_text(size = 8.2),
    axis.text = element_text(size = 7.1, colour = "black"),
    plot.margin = margin(2, 10, 2, 4)
  )

ggsave(file.path(OUT, "Fig4d.pdf"), pD, width = 7.2, height = 4.2, units = "in", device = "pdf")
message("Wrote: ", file.path(OUT, "Fig4d.pdf"))

writeLines(capture.output(sessionInfo()), "../environment/R_sessionInfo_06_pam_light_response.txt")
