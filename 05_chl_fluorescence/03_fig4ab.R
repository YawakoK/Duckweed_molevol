# Step 03: Fig. 4a (Fv/Fm time course, mean +/- SE) and Fig. 4b (Day 1 vs Day 4 within
#          lineage, brackets = Tukey-adjusted contrasts from Supplementary Data 12).
#
# The panels are drawn at the size they occupy in the combined Fig. 4 (7.2 in wide; top row
# split 0.85 : 1.15).  The combined Fig. 4 places these two panels above panels c-d from
# ../06_pam_light_response; the final layout (species abbreviations, font weights) was
# finished in a vector editor.
#
# Input : output/SupplementaryData10_FvFm_extracted_data.csv                      (step 01)
#         output/SupplementaryData12_FvFm_time_pairwise_comparisons_within_species.csv (step 02)
# Output: output/Fig4a.pdf, output/Fig4b.pdf
# Run from this directory:  Rscript 03_fig4ab.R
# Author: Natsu Katayama

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
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

fig4_samples <- c("Sp7498", "Lp9387", "Lgp8L", "LaM10E", "Wh9525", "Wa")
fig4_sample_labels <- c(
  "Sp7498" = "Spirodela polyrhiza",
  "Lp9387" = "Landoltia punctata",
  "Lgp8L" = "Lemna gibba",
  "LaM10E" = "Lemna aequinoctialis",
  "Wh9525" = "Wolffiella hyalina",
  "Wa" = "Wolffia australiana"
)
fig4_species_cols <- c(
  "Sp7498" = "#8c77af",
  "Lp9387" = "#6f5a9a",
  "Lgp8L" = "#b39ddb",
  "LaM10E" = "#9bb85c",
  "Wh9525" = "#4f7f3a",
  "Wa" = "#73a244"
)
fig4_species_shapes <- c(
  "Sp7498" = 16,
  "Lp9387" = 17,
  "Lgp8L" = 15,
  "LaM10E" = 18,
  "Wh9525" = 3,
  "Wa" = 8
)

fvfm <- read_csv(file.path(OUT, "SupplementaryData10_FvFm_extracted_data.csv"), show_col_types = FALSE) |>
  mutate(
    sample = factor(sample, levels = fig4_samples),
    time = factor(time, levels = paste0("DAY", 0:4))
  )

fvfm_summary <- fvfm |>
  filter(time %in% c("DAY0", "DAY1", "DAY4")) |>
  group_by(sample, time) |>
  summarise(
    mean_FvFm = mean(Fv_Fm, na.rm = TRUE),
    se_FvFm = sd(Fv_Fm, na.rm = TRUE) / sqrt(sum(!is.na(Fv_Fm))),
    .groups = "drop"
  )

fvfm_analysis <- fvfm |>
  filter(time %in% c("DAY1", "DAY4"))

fvfm_recovery_stats <- read_csv(
  file.path(OUT, "SupplementaryData12_FvFm_time_pairwise_comparisons_within_species.csv"),
  show_col_types = FALSE
) |>
  filter(contrast == "DAY1 - DAY4", sample %in% fig4_samples) |>
  mutate(
    sample = factor(sample, levels = fig4_samples),
    p_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "n.s."
    )
  ) |>
  left_join(
    fvfm_summary |>
      filter(time %in% c("DAY1", "DAY4")) |>
      group_by(sample) |>
      summarise(y = max(mean_FvFm + se_FvFm, na.rm = TRUE) + 0.035, .groups = "drop"),
    by = "sample"
  ) |>
  mutate(
    x = as.numeric(sample),
    x_start = x - 0.18,
    x_end = x + 0.18,
    y_tip = y - 0.012,
    y_text = y + 0.007
  )

pA <- ggplot(fvfm_summary, aes(x = time, y = mean_FvFm, color = sample, shape = sample, group = sample)) +
  geom_line(linewidth = 0.50) +
  geom_point(size = 1.75, stroke = 0.55) +
  geom_errorbar(
    aes(ymin = mean_FvFm - se_FvFm, ymax = mean_FvFm + se_FvFm),
    width = 0.07,
    linewidth = 0.30
  ) +
  scale_color_manual(
    values = fig4_species_cols,
    breaks = fig4_samples,
    labels = fig4_sample_labels[fig4_samples],
    name = "Species"
  ) +
  scale_shape_manual(
    values = fig4_species_shapes,
    breaks = fig4_samples,
    labels = fig4_sample_labels[fig4_samples],
    name = "Species"
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 0, size = 1.8)),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_x_discrete(expand = expansion(add = c(0.10, 0.22))) +
  coord_cartesian(ylim = c(0.45, 0.85), clip = "off") +
  labs(x = "Day", y = expression(F[v]/F[m])) +
  theme_base +
  theme(
    axis.title = element_text(size = 7.8),
    axis.text = element_text(size = 6.9, colour = "black"),
    legend.position = "top",
    legend.text = element_text(size = 5.4, face = "italic"),
    legend.title = element_text(size = 6.2),
    legend.key.width = unit(3.0, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 4, 2, 4)
  )

set.seed(1)   # position_jitterdodge() for the individual points in Fig. 4b
pB <- ggplot(fvfm_analysis, aes(x = sample, y = Fv_Fm)) +
  geom_boxplot(
    aes(fill = time),
    width = 0.48,
    outlier.shape = NA,
    colour = "black",
    linewidth = 0.35,
    position = position_dodge(width = 0.58)
  ) +
  geom_point(
    aes(fill = time),
    shape = 21,
    size = 1.45,
    alpha = 0.65,
    colour = "black",
    stroke = 0.15,
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.58)
  ) +
  geom_segment(
    data = fvfm_recovery_stats,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  geom_segment(
    data = fvfm_recovery_stats,
    aes(x = x_start, xend = x_start, y = y, yend = y_tip),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  geom_segment(
    data = fvfm_recovery_stats,
    aes(x = x_end, xend = x_end, y = y, yend = y_tip),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  geom_text(
    data = fvfm_recovery_stats,
    aes(x = x, y = y_text, label = p_label),
    inherit.aes = FALSE,
    size = 2.25
  ) +
  scale_x_discrete(
    labels = c(
      "Sp7498" = "Spirodela\npolyrhiza",
      "Lp9387" = "Landoltia\npunctata",
      "Lgp8L" = "Lemna\ngibba",
      "LaM10E" = "Lemna\naequinoctialis",
      "Wh9525" = "Wolffiella\nhyalina",
      "Wa" = "Wolffia\naustraliana"
    )
  ) +
  scale_fill_manual(values = c("DAY1" = "#D8E6F3", "DAY4" = "#4F8BB8"),
                    labels = c("Day 1", "Day 4"), name = NULL) +
  coord_cartesian(ylim = c(0.45, 0.78), clip = "off") +
  labs(x = "Species", y = expression(F[v]/F[m])) +
  theme_base +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 7.8),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = "italic", size = 6.1),
    axis.text.y = element_text(size = 6.8, colour = "black"),
    legend.text = element_text(size = 6.6),
    plot.margin = margin(2, 6, 2, 4)
  )

# top row of the combined Fig. 4: 7.2 in x (0.88 / 3.32 of 10.8 in)
row_h <- 10.8 * 0.88 / (0.88 + 1.10 + 1.34)
ggsave(file.path(OUT, "Fig4a.pdf"), pA, width = 7.2 * 0.85 / 2, height = row_h, units = "in")
ggsave(file.path(OUT, "Fig4b.pdf"), pB, width = 7.2 * 1.15 / 2, height = row_h, units = "in")

cat("Fig. 4b brackets (Day 1 vs Day 4, Tukey):\n")
print(fvfm_recovery_stats |> select(sample, p.value, p_label), row.names = FALSE)
# -> Fig. 4b: Sp ***, Lp n.s., Lg *, La n.s., Wh n.s., Wa **
