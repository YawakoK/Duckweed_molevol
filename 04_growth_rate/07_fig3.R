# Step 07: Fig. 3 (a: RGR by species and light; b: RGR by anthocyanin status and light).
#
# Boxplots of the 72 well-level Day 0 - Day 4 RGR values; annotations are the GLMM results
# of step 06 (a: HL/LL emmeans ratio and P within species, species x light P;
# b: non-antho/antho emmeans ratio and P within light, anthocyanin status x light P).
# Point jitter uses set.seed(1); the submitted figure was additionally edited in a vector
# editor (species abbreviations, italics, layout), data and annotations unchanged.
#
# Input : output/SupplementaryData06_RGR_by_well_Day0_Day4.csv (04)
#         output/Fig3a_HL_LL_ratio_by_species.csv, output/Fig3a_species_light_TypeII_Wald.csv (06)
#         output/SupplementaryData07_RGR_anthocyanin_light_TypeII_Wald.csv (06)
#         output/SupplementaryData09_RGR_pairwise_anthocyanin_contrasts_by_light.csv (06)
# Output: output/Fig3.pdf
#         ../environment/R_sessionInfo_04_growth_rate.txt
# Run from this directory:  Rscript 07_fig3.R
# Author: Natsu Katayama
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

set.seed(1)

fig3_samples <- c("Sp", "Lp", "Lgp8L", "M10E", "Wh", "Wa")
light_levels <- c("PFD100", "PFD1000")
light_cols <- c("PFD100" = "#00BFC4", "PFD1000" = "#F8766D")
light_labels <- c("PFD100" = "LL", "PFD1000" = "HL")
species_cols <- c(
  "Sp" = "#8c77af",
  "Lp" = "#6f5a9a",
  "Lgp8L" = "#b39ddb",
  "M10E" = "#9bb85c",
  "Wh" = "#4f7f3a",
  "Wa" = "#73a244"
)
species_shapes <- c(
  "Sp" = 16,
  "Lp" = 17,
  "Lgp8L" = 15,
  "M10E" = 18,
  "Wh" = 3,
  "Wa" = 8
)
species_axis_labels <- c(
  "Sp" = "Spirodela\npolyrhiza",
  "Lp" = "Landoltia\npunctata",
  "Lgp8L" = "Lemna\ngibba",
  "M10E" = "Lemna\naequinoctialis",
  "Wh" = "Wolffiella\nhyalina",
  "Wa" = "Wolffia\naustraliana"
)
species_legend_labels <- c(
  "Sp" = "Spirodela\npolyrhiza",
  "Lp" = "Landoltia\npunctata",
  "Lgp8L" = "Lemna\ngibba",
  "M10E" = "Lemna\naequinoctialis",
  "Wh" = "Wolffiella\nhyalina",
  "Wa" = "Wolffia\naustraliana"
)
antho_cols <- c("antho" = "#8c77af", "non-antho" = "#73a244")
antho_labels <- c(
  "antho" = "Accumulating",
  "non-antho" = "Non-accumulating"
)

p_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 1e-4 ~ "****",
    p < 1e-3 ~ "***",
    p < 1e-2 ~ "**",
    p < 5e-2 ~ "*",
    TRUE ~ "n.s."
  )
}

p_text <- function(p) {
  ifelse(p < 1e-4, "p < 1e-4", paste0("p = ", signif(p, 2)))
}

rgr <- read_csv("output/SupplementaryData06_RGR_by_well_Day0_Day4.csv", show_col_types = FALSE) |>
  mutate(
    sample = factor(sample, levels = fig3_samples),
    light = factor(light, levels = light_levels),
    antho = factor(antho, levels = c("antho", "non-antho"))
  )

species_y <- rgr |>
  group_by(sample) |>
  summarise(y = max(rgr_per_day, na.rm = TRUE) + 0.045, .groups = "drop")

light_contrasts <- read_csv("output/Fig3a_HL_LL_ratio_by_species.csv", show_col_types = FALSE) |>
  mutate(
    sample = factor(sample, levels = fig3_samples),
    x = as.numeric(sample),
    label = p_label(p.value)
  ) |>
  left_join(species_y, by = "sample")

sample_interaction_p <- read_csv("output/Fig3a_species_light_TypeII_Wald.csv", show_col_types = FALSE) |>
  filter(Effect == "sample:light") |>
  pull(`Pr(>Chisq)`)

antho_y <- rgr |>
  group_by(light) |>
  summarise(y = max(rgr_per_day, na.rm = TRUE) * 1.12, .groups = "drop")

pairwise_antho <- read_csv("output/SupplementaryData09_RGR_pairwise_anthocyanin_contrasts_by_light.csv",
                           show_col_types = FALSE) |>
  mutate(
    light = factor(light, levels = light_levels),
    x = as.numeric(light),
    label = p_label(p.value)
  ) |>
  left_join(antho_y, by = "light")

antho_interaction_p <- read_csv("output/SupplementaryData07_RGR_anthocyanin_light_TypeII_Wald.csv",
                                show_col_types = FALSE) |>
  filter(Effect == "antho:light") |>
  pull(`Pr(>Chisq)`)

dodge_width <- 0.75

theme_fig3 <- theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 8),
    legend.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(4, "mm"),
    plot.tag = element_text(size = 11, face = "bold"),
    plot.tag.position = c(0.01, 0.98),
    plot.margin = margin(3, 4, 3, 4)
  )

species_anno <- light_contrasts |>
  mutate(
    x1 = x - dodge_width / 4,
    x2 = x + dodge_width / 4,
    hl_vs_ll_ratio = 1 / ratio,
    ratio_label = paste0(label, "\nratio = ", sprintf("%.2f", hl_vs_ll_ratio)),
    y_text = y + 0.030
  )

p_species <- ggplot(rgr, aes(x = sample, y = rgr_per_day, fill = light)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.68,
    alpha = 0.72,
    colour = "black",
    linewidth = 0.35,
    position = position_dodge(width = dodge_width)
  ) +
  geom_point(
    colour = "grey35",
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = dodge_width),
    size = 1.1,
    alpha = 0.72
  ) +
  geom_segment(
    data = species_anno,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  geom_segment(
    data = species_anno,
    aes(x = x1, xend = x1, y = y * 0.985, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  geom_segment(
    data = species_anno,
    aes(x = x2, xend = x2, y = y * 0.985, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  geom_text(
    data = species_anno,
    aes(x = x, y = y_text, label = ratio_label),
    inherit.aes = FALSE,
    size = 2.0,
    lineheight = 0.82
  ) +
  annotate(
    "text",
    x = 3.5,
    y = max(species_anno$y_text, na.rm = TRUE) + 0.065,
    label = paste0("species x light: ", p_text(sample_interaction_p)),
    size = 2.2
  ) +
  scale_fill_manual(values = light_cols, labels = light_labels, name = NULL) +
  scale_x_discrete(labels = species_axis_labels) +
  labs(x = NULL, y = expression(RGR~(day^-1))) +
  coord_cartesian(ylim = c(0.05, max(species_anno$y_text, na.rm = TRUE) + 0.115), clip = "off") +
  theme_fig3 +
  theme(
    axis.text.x = element_text(face = "italic", lineheight = 0.85, size = 7),
    legend.margin = margin(0, 0, 0, 0)
  )

antho_anno <- pairwise_antho |>
  mutate(
    x1 = x - dodge_width / 4,
    x2 = x + dodge_width / 4,
    y_text = y * 1.08,
    ratio_label = paste0(label, "\nratio = ", sprintf("%.2f", ratio))
  )

p_antho <- ggplot(rgr, aes(x = light, y = rgr_per_day, fill = antho)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65,
    alpha = 0.72,
    colour = "black",
    linewidth = 0.35,
    position = position_dodge(width = dodge_width)
  ) +
  geom_point(
    aes(colour = sample, shape = sample, group = antho),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = dodge_width),
    size = 1.35,
    alpha = 0.82
  ) +
  geom_segment(
    data = antho_anno,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  geom_segment(
    data = antho_anno,
    aes(x = x1, xend = x1, y = y * 0.985, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  geom_segment(
    data = antho_anno,
    aes(x = x2, xend = x2, y = y * 0.985, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  geom_text(
    data = antho_anno,
    aes(x = x, y = y_text, label = ratio_label),
    inherit.aes = FALSE,
    size = 2.1,
    lineheight = 0.82
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(antho_anno$y_text, na.rm = TRUE) * 1.10,
    label = paste0("anthocyanin status x light: ", p_text(antho_interaction_p)),
    size = 2.2
  ) +
  scale_fill_manual(values = antho_cols, labels = antho_labels, name = NULL, guide = "none") +
  scale_colour_manual(values = species_cols, labels = species_legend_labels, name = NULL) +
  scale_shape_manual(values = species_shapes, labels = species_legend_labels, name = NULL) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_x_discrete(labels = light_labels) +
  labs(x = "Light condition", y = expression(RGR~(day^-1))) +
  coord_cartesian(ylim = c(0.05, max(antho_anno$y_text, na.rm = TRUE) * 1.16), clip = "off") +
  theme_fig3 +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 5.2, lineheight = 0.82, face = "italic"),
    legend.margin = margin(0, 0, 0, 0)
  )

fig3 <- (p_species + p_antho + plot_layout(widths = c(1.45, 1))) +
  plot_annotation(tag_levels = "A") &
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("output/Fig3.pdf", fig3, width = 7.2, height = 3.7, units = "in", device = "pdf")

writeLines(capture.output(sessionInfo()), "../environment/R_sessionInfo_04_growth_rate.txt")
