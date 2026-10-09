# Step 02: PCA overview of the three metabolome datasets (Supplementary Fig. 8).
#
# Input : output/metabolite_abundance_long.csv (step 01)
# Output: output/SupplementaryFig8.pdf
#         output/SupplementaryFig8_PCA_scores.csv (PC1/PC2 scores and % variance)
# Run from this directory: Rscript 02_pca.R
# Author: Natsu Katayama
#
# For PCA only, abundances are transformed as log2(value + 1); every metabolite
# is centred and scaled to unit variance (prcomp, center = TRUE, scale. = TRUE).

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

abundance_long <- read_csv("output/metabolite_abundance_long.csv", show_col_types = FALSE)

run_pca <- function(platform) {
  mat <- abundance_long %>%
    filter(Platform == platform) %>%
    select(Metabolite, sample, value) %>%
    pivot_wider(names_from = sample, values_from = value) %>%
    tibble::column_to_rownames("Metabolite") %>%
    as.matrix()
  mat <- t(log2(mat + 1))                      # rows = samples, columns = metabolites
  prcomp(mat, center = TRUE, scale. = TRUE)
}

pca_plot <- function(res, title) {
  scores <- as.data.frame(res$x) %>%
    tibble::rownames_to_column("sample") %>%
    separate(sample, into = c("species", "light", "rep"), sep = "_", remove = FALSE)
  ve <- summary(res)$importance[2, ]
  ggplot(scores, aes(PC1, PC2, color = light, shape = species)) +
    geom_point(size = 1.5) +
    theme_bw() +
    labs(
      title = title,
      x = paste0("PC1 (", round(100 * ve[1], 1), "%)"),
      y = paste0("PC2 (", round(100 * ve[2], 1), "%)")
    )
}

pca_lcms <- run_pca("LC-MS/MS")
pca_gcms <- run_pca("GC-MS/MS")
pca_cems <- run_pca("CE-TOF MS")

p_all <- pca_plot(pca_lcms, "LC-MS PCA") |
  pca_plot(pca_gcms, "GC-MS PCA") |
  pca_plot(pca_cems, "CE-MS PCA")

ggsave("output/SupplementaryFig8.pdf", p_all, width = 210, height = 60, units = "mm")

scores_out <- bind_rows(lapply(
  list("LC-MS/MS" = pca_lcms, "GC-MS/MS" = pca_gcms, "CE-TOF MS" = pca_cems),
  function(res) {
    ve <- summary(res)$importance[2, ]
    tibble(
      sample = rownames(res$x), PC1 = res$x[, 1], PC2 = res$x[, 2],
      PC1_percent_variance = 100 * ve[1], PC2_percent_variance = 100 * ve[2]
    )
  }
), .id = "Platform")
write_csv(scores_out, "output/SupplementaryFig8_PCA_scores.csv")

# -> Supplementary Fig. 8 axis labels: % variance explained by PC1 / PC2
print(scores_out %>% distinct(Platform, PC1_percent_variance, PC2_percent_variance))
