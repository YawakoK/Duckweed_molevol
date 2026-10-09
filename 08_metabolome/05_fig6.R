# Step 05: Fig. 6a-e - category log2FC scatter plots (S. polyrhiza vs
#          W. australiana) and log2 abundance of three representative
#          metabolites per category.
#
# Input : output/SupplementaryData21_all_metabolite_log2FC_statistics.csv (step 03)
#         output/SupplementaryData22_category_metabolite_log2FC_statistics.csv (step 04)
#         output/metabolite_log2_abundance_long.csv (step 03)
# Output: output/Fig6.pdf (R panels; the submitted Fig. 6 is the same panels
#         re-arranged in Affinity Designer)
#         output/Fig6_representative_metabolites_source_data.csv
# Run from this directory: Rscript 05_fig6.R
# Author: Natsu Katayama
#
# Scatter colours: purple = Sp > Wa, green = Wa > Sp (interaction FDR < 0.05),
# grey = not significant. Abundance panels: replicate points, lineage means
# joined by lines; the label is the interaction FDR (Supplementary Data 21).

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(ggrepel)

supp21 <- read_csv("output/SupplementaryData21_all_metabolite_log2FC_statistics.csv",
                   show_col_types = FALSE)
supp22 <- read_csv("output/SupplementaryData22_category_metabolite_log2FC_statistics.csv",
                   show_col_types = FALSE)
log2_long <- read_csv("output/metabolite_log2_abundance_long.csv", show_col_types = FALSE)

# Panel definitions (row order of the figure)
category_map <- tibble(
  panel_assignment = c(
    "Fig. 6A: Flavonoids",
    "Fig. 6B: Amino acids / peptides",
    "Fig. 6C: Central carbon metabolism",
    "Fig. 6D: Carbohydrates / sugar metabolism",
    "Fig. 6E: Energy / redox"
  ),
  CategoryTitle = c(
    "A Flavonoids",
    "B Amino acids / peptides",
    "C Central carbon metabolism",
    "D Carbohydrates / sugar metabolism",
    "E Energy / redox"
  )
)

# Three representative metabolites per category, chosen for biological
# interpretability among the significant species-biased metabolites.
# SelectedLabel = label used in the scatter plot.
representatives <- tribble(
  ~panel_assignment,                              ~Platform,   ~Metabolite,                                 ~SelectedLabel,
  "Fig. 6A: Flavonoids",                          "LC-MS/MS",  "Nevadensin 7-rutinoside",                   "Nevadensin 7-rutinoside",
  "Fig. 6A: Flavonoids",                          "LC-MS/MS",  "Apigenin 7-(6''-malonylneohesperidoside)",  "Apigenin glycoside",
  "Fig. 6A: Flavonoids",                          "LC-MS/MS",  "Flavonol 7-O-beta-D-glucoside;",            "Flavonol glucoside",
  "Fig. 6B: Amino acids / peptides",              "GC-MS/MS",  "Glycine-3TMS",                              "Glycine",
  "Fig. 6B: Amino acids / peptides",              "GC-MS/MS",  "Serine-3TMS",                               "Serine",
  "Fig. 6B: Amino acids / peptides",              "GC-MS/MS",  "N-Acetylglutamine-3TMS",                    "N-Acetylglutamine",
  "Fig. 6C: Central carbon metabolism",           "GC-MS/MS",  "Citric acid-4TMS",                          "Citric acid",
  "Fig. 6C: Central carbon metabolism",           "GC-MS/MS",  "Aconitic acid-3TMS",                        "Aconitic acid",
  "Fig. 6C: Central carbon metabolism",           "GC-MS/MS",  "Malic acid-3TMS",                           "Malic acid",
  "Fig. 6D: Carbohydrates / sugar metabolism",    "GC-MS/MS",  "Glucose-meto-5TMS(2)",                      "Glucose",
  "Fig. 6D: Carbohydrates / sugar metabolism",    "GC-MS/MS",  "Trehalose-8TMS",                            "Trehalose",
  "Fig. 6D: Carbohydrates / sugar metabolism",    "GC-MS/MS",  "Sucrose-8TMS",                              "Sucrose",
  "Fig. 6E: Energy / redox",                      "CE-TOF MS", "T-049+_ATP[M-2H]",                          "ATP",
  "Fig. 6E: Energy / redox",                      "CE-TOF MS", "T-048+_ADP[M-2H]",                          "ADP",
  "Fig. 6E: Energy / redox",                      "CE-TOF MS", "T-060_NAD",                                 "NAD"
) %>%
  mutate(order = row_number()) %>%
  left_join(supp22 %>% select(MS, Metabolite, MetLabel), by = c("Platform" = "MS", "Metabolite")) %>%
  left_join(supp21 %>% select(Platform, Metabolite, log2FC_Sp, log2FC_Wa, p_interaction, FDR, sig_class),
            by = c("Platform", "Metabolite"))
stopifnot(!any(is.na(representatives$MetLabel)), !any(is.na(representatives$FDR)))

# -> Results: interaction FDR of the representative metabolites (labels in Fig. 6a-e)
print(representatives %>% select(MetLabel, log2FC_Sp, log2FC_Wa, FDR, sig_class), n = Inf)

abundance_df <- log2_long %>%
  inner_join(representatives %>% select(panel_assignment, order, Platform, Metabolite, MetLabel),
             by = c("Platform", "Metabolite")) %>%
  mutate(
    Lineage = factor(Lineage, levels = c("Sp", "Wa")),
    Light = factor(Light, levels = c("LL", "HL")),
    MetLabel = factor(MetLabel, levels = representatives$MetLabel)
  ) %>%
  arrange(order, Lineage, Light, sample)

write_csv(
  abundance_df %>%
    left_join(representatives %>% select(Platform, Metabolite, log2FC_Sp, log2FC_Wa, p_interaction, FDR, sig_class),
              by = c("Platform", "Metabolite")) %>%
    select(panel_assignment, Platform, Metabolite, MetLabel, sample, rep, Lineage, Light,
           value, pseudo, log2y, log2FC_Sp, log2FC_Wa, p_interaction, FDR, sig_class),
  "output/Fig6_representative_metabolites_source_data.csv"
)

scatter_df <- supp22 %>%
  left_join(category_map, by = "panel_assignment") %>%
  left_join(representatives %>% select(MS = Platform, Metabolite, SelectedLabel),
            by = c("MS", "Metabolite")) %>%
  mutate(
    CategoryTitle = factor(CategoryTitle, levels = category_map$CategoryTitle),
    plot_class = factor(plot_class, levels = c("Sp > Wa", "Wa > Sp", "Not significant / other")),
    is_selected = !is.na(SelectedLabel)
  )

# ---- plotting (style as submitted) -------------------------------------------------
plot_class_cols <- c(
  "Sp > Wa" = "#8c77af",
  "Wa > Sp" = "#73a244",
  "Not significant / other" = "grey70"
)
lineage_cols <- c("Sp" = "#8c77af", "Wa" = "#73a244")
lineage_labs <- c("Sp" = "Spirodela polyrhiza", "Wa" = "Wolffia australiana")

format_fdr <- function(x) {
  ifelse(is.na(x), "FDR = NA", ifelse(x < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", x)))
}

# common axis range over all five scatter panels
scatter_limits <- range(c(scatter_df$Sp_HL_response, scatter_df$Wa_HL_response), na.rm = TRUE)
scatter_pad <- diff(scatter_limits) * 0.06
scatter_limits <- scatter_limits + c(-scatter_pad, scatter_pad)

base_theme <- theme_bw(base_size = 6.8) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.18, colour = "grey92"),
    axis.text = element_text(colour = "black", size = 5.6),
    axis.title = element_text(size = 6.3),
    strip.background = element_blank(),
    strip.text = element_text(size = 5.8, lineheight = 0.88),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 6.9, hjust = 0, lineheight = 0.88,
                              margin = margin(b = -4.0)),
    plot.margin = margin(-2.0, 1.5, 8.0, 1.5)
  )

make_scatter <- function(panel, title) {
  dat <- scatter_df %>% filter(panel_assignment == panel)
  ggplot(dat, aes(x = Sp_HL_response, y = Wa_HL_response)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey75") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.28, linetype = "dashed", colour = "grey55") +
    geom_point(aes(colour = plot_class), size = 1.05, alpha = 0.85) +
    geom_point(data = dat %>% filter(is_selected), shape = 21, fill = NA,
               colour = "black", stroke = 0.35, size = 1.85) +
    geom_text_repel(
      data = dat %>% filter(is_selected),
      aes(label = SelectedLabel),
      size = 1.65, lineheight = 0.82, colour = "grey15",
      box.padding = 0.12, point.padding = 0.10, min.segment.length = 0,
      segment.colour = "grey35", segment.size = 0.15,
      max.overlaps = Inf, seed = 7
    ) +
    coord_equal(xlim = scatter_limits, ylim = scatter_limits) +
    scale_colour_manual(values = plot_class_cols, drop = FALSE) +
    labs(title = title, x = "log2FC in S. polyrhiza", y = "log2FC in W. australiana") +
    base_theme +
    theme(legend.position = "none")
}

make_abundance <- function(panel) {
  stats_block <- representatives %>% filter(panel_assignment == panel)
  met_levels <- stats_block$MetLabel
  dat <- abundance_df %>%
    filter(panel_assignment == panel) %>%
    mutate(MetLabel = factor(as.character(MetLabel), levels = met_levels))
  means <- dat %>%
    group_by(MetLabel, Lineage, Light) %>%
    summarise(mean_log2 = mean(log2y, na.rm = TRUE), .groups = "drop")
  pvals <- stats_block %>%
    mutate(
      MetLabel = factor(MetLabel, levels = met_levels),
      Light = factor("LL", levels = c("LL", "HL")),
      label = format_fdr(FDR)
    )

  ggplot(means, aes(x = Light, y = mean_log2, group = Lineage, colour = Lineage)) +
    geom_point(data = dat, aes(x = Light, y = log2y, colour = Lineage), inherit.aes = FALSE,
               position = position_jitter(width = 0.08, height = 0, seed = 1),
               size = 0.45, alpha = 0.35) +
    geom_line(linewidth = 0.35, alpha = 0.9) +
    geom_point(size = 0.95, alpha = 0.95) +
    geom_text(data = pvals, aes(x = Light, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.15, vjust = 1.4, size = 1.45, colour = "grey20") +
    facet_wrap(~MetLabel, scales = "free_y", ncol = min(3, length(met_levels)),
               labeller = label_wrap_gen(width = 18)) +
    scale_colour_manual(values = lineage_cols, labels = lineage_labs, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    labs(x = NULL, y = "log2 abundance", colour = NULL) +
    base_theme +
    theme(legend.position = "none", axis.title.y = element_text(size = 6.0),
          panel.spacing = unit(0.35, "lines"))
}

rows <- lapply(seq_len(nrow(category_map)), function(i) {
  make_scatter(category_map$panel_assignment[[i]], category_map$CategoryTitle[[i]]) +
    make_abundance(category_map$panel_assignment[[i]]) +
    plot_layout(widths = c(0.92, 2.35))
})

legend_df <- tibble(
  x = c(0.10, 0.28, 0.42, 0.63, 0.82),
  y = 0.5,
  label = c("Not significant / other", "Sp > Wa", "Wa > Sp",
            "Spirodela polyrhiza", "Wolffia australiana"),
  type = c("scatter", "scatter", "scatter", "lineage", "lineage"),
  colour = c("grey70", "#8c77af", "#73a244", "#8c77af", "#73a244")
)

legend_plot <- ggplot() +
  annotate("text", x = 0.00, y = 0.5, label = "log2FC class", hjust = 0, size = 1.95) +
  geom_point(data = legend_df %>% filter(type == "scatter"),
             aes(x = x, y = y, colour = label), size = 1.7) +
  annotate("text", x = 0.54, y = 0.5, label = "Lineage", hjust = 0, size = 1.95) +
  geom_segment(data = legend_df %>% filter(type == "lineage"),
               aes(x = x - 0.014, xend = x + 0.014, y = y, yend = y, colour = label),
               linewidth = 0.38) +
  geom_point(data = legend_df %>% filter(type == "lineage"),
             aes(x = x, y = y, colour = label), size = 1.5) +
  geom_text(data = legend_df, aes(x = x + 0.025, y = y, label = label), hjust = 0, size = 1.75) +
  scale_colour_manual(values = setNames(legend_df$colour, legend_df$label), guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void()

combined <- wrap_plots(
  plotlist = list(wrap_plots(rows, ncol = 1), legend_plot),
  ncol = 1,
  heights = c(1, 0.032)
) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("output/Fig6.pdf", combined, width = 7.2, height = 11.4, units = "in", device = "pdf")

writeLines(capture.output(sessionInfo()), "../environment/R_sessionInfo_08_metabolome.txt")
