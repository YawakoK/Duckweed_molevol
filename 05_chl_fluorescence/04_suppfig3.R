# Step 04: Supplementary Fig. 3 - Fv/Fm distributions of the six lineages at Day 0, 1 and 4,
#          with compact letters among lineages within each day (Supplementary Data 14,
#          Sidak-adjusted). Box fill = anthocyanin status; points = individual measurements.
#
# Input : output/SupplementaryData10_FvFm_extracted_data.csv                   (step 01)
#         output/SupplementaryData11_FvFm_summary_by_species_and_time.csv      (step 02)
#         output/SupplementaryData14_FvFm_species_group_letters_within_time.csv (step 02)
# Output: output/SupplementaryFig3.pdf
#         ../environment/R_sessionInfo_05_chl_fluorescence.txt
# Run from this directory:  Rscript 04_suppfig3.R
# Author: Natsu Katayama

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

OUT <- "output"
analysis_samples <- c("Sp7498", "Lp9387", "Lgp8L", "LaM10E", "Wh9525", "Wa")
analysis_times <- c("DAY0", "DAY1", "DAY4")

axis_labels <- c(
  "Sp7498" = "Spirodela\npolyrhiza",
  "Lp9387" = "Landoltia\npunctata",
  "Lgp8L" = "Lemna\ngibba",
  "LaM10E" = "Lemna\naequinoctialis",
  "Wh9525" = "Wolffiella\nhyalina",
  "Wa" = "Wolffia\naustraliana"
)
species_cols <- c(
  "Sp7498" = "#1F4E99",
  "Lp9387" = "#2F7FC1",
  "Lgp8L" = "#69BDE5",
  "LaM10E" = "#F2C94C",
  "Wh9525" = "#E07B24",
  "Wa" = "#D95F02"
)
antho_cols <- c("antho" = "#8c77af", "non_antho" = "#73a244")
antho_labels <- c("antho" = "Anthocyanin-present", "non_antho" = "Anthocyanin-absent")

theme_fvfm <- theme_bw(base_size = 10) +
  theme(
    panel.border = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "top"
  )
species_axis_theme <- theme(
  axis.text.x = element_text(angle = 35, hjust = 0.65, vjust = 1, face = "italic",
                             colour = "black", size = 8.5, margin = margin(t = 20)),
  plot.margin = margin(5.5, 8, 54, 22)
)

fac <- function(d) {
  d |> mutate(sample = factor(sample, levels = analysis_samples),
              time = factor(time, levels = analysis_times))
}
fvfm <- read_csv(file.path(OUT, "SupplementaryData10_FvFm_extracted_data.csv"), show_col_types = FALSE) |>
  filter(time %in% analysis_times) |>
  fac() |>
  mutate(antho = factor(antho, levels = c("antho", "non_antho")))
summary_df <- read_csv(file.path(OUT, "SupplementaryData11_FvFm_summary_by_species_and_time.csv"),
                       show_col_types = FALSE) |> fac()
letter_plot_df <- read_csv(file.path(OUT, "SupplementaryData14_FvFm_species_group_letters_within_time.csv"),
                           show_col_types = FALSE) |>
  fac() |>
  select(time, sample, .group) |>
  left_join(summary_df, by = c("time", "sample")) |>
  mutate(label_y = pmin(mean_FvFm + se_FvFm + 0.035, 0.88))

set.seed(1)   # geom_jitter()
p_box <- ggplot(fvfm, aes(x = sample, y = Fv_Fm, fill = antho)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, colour = "black") +
  geom_jitter(aes(color = sample), width = 0.12, size = 2.2, alpha = 0.45) +
  geom_text(
    data = letter_plot_df,
    aes(x = sample, y = label_y, label = .group),
    inherit.aes = FALSE,
    size = 4.2
  ) +
  facet_wrap(~ time, nrow = 1) +
  scale_fill_manual(values = antho_cols, labels = antho_labels, name = NULL) +
  scale_color_manual(values = species_cols, guide = "none") +
  scale_x_discrete(labels = axis_labels[analysis_samples]) +
  coord_cartesian(ylim = c(0.45, 0.85)) +
  labs(x = "Species", y = expression(F[v]/F[m])) +
  theme_fvfm +
  species_axis_theme

ggsave(file.path(OUT, "SupplementaryFig3.pdf"), p_box, width = 10.5, height = 6.0)

# record the versions of the step-02 statistics packages as well
invisible(lapply(c("car", "emmeans", "multcomp"), loadNamespace))
writeLines(capture.output(sessionInfo()),"../environment/R_sessionInfo_05_chl_fluorescence.txt")
