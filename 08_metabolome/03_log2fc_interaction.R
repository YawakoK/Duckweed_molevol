# Step 03: Pseudocount, log2 transformation, within-species log2FC and the
#          Lineage x Light interaction test for every metabolite
#          (Supplementary Data 21).
#
# Input : output/metabolite_abundance_long.csv (step 01)
# Output: output/metabolite_log2_abundance_long.csv (adds pseudo, log2y)
#         output/SupplementaryData21_all_metabolite_log2FC_statistics.csv
# Run from this directory: Rscript 03_log2fc_interaction.R
# Author: Natsu Katayama
#
# Methods ("Metabolome analysis"):
#   pseudo = half of the smallest positive value of each metabolite (1e-9 / 2
#            if a metabolite has no positive value); log2y = log2(value + pseudo)
#   log2FC (per species) = mean log2y under HL - mean log2y under LL
#   lm(log2y ~ Lineage * Light); p_interaction = anova() F test of Lineage:Light
#   FDR = Benjamini-Hochberg across metabolites, separately within each platform
#   sig_class: FDR < 0.05 and log2FC_Wa > log2FC_Sp -> "Wa > Sp";
#              FDR < 0.05 and log2FC_Sp > log2FC_Wa -> "Sp > Wa"

library(dplyr)
library(tidyr)
library(readr)

abundance_long <- read_csv("output/metabolite_abundance_long.csv", show_col_types = FALSE) %>%
  mutate(
    Lineage = factor(Lineage, levels = c("Sp", "Wa")),
    Light = factor(Light, levels = c("LL", "HL"))
  )

log2_long <- abundance_long %>%
  group_by(Platform, Metabolite) %>%
  mutate(
    min_pos = suppressWarnings(min(value[value > 0], na.rm = TRUE)),
    min_pos = ifelse(is.finite(min_pos) & min_pos > 0, min_pos, 1e-9),
    pseudo = min_pos / 2,
    log2y = log2(value + pseudo)
  ) %>%
  ungroup() %>%
  select(-min_pos)
write_csv(log2_long, "output/metabolite_log2_abundance_long.csv")

log2fc <- log2_long %>%
  group_by(Platform, Ion, Metabolite, Lineage, Light) %>%
  summarise(mean_log2y = mean(log2y, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = c(Light, Lineage),
    values_from = mean_log2y,
    names_glue = "log2{Light}_mean_{Lineage}"
  ) %>%
  mutate(
    log2FC_Sp = log2HL_mean_Sp - log2LL_mean_Sp,
    log2FC_Wa = log2HL_mean_Wa - log2LL_mean_Wa
  )

interaction_p <- log2_long %>%
  group_by(Platform, Metabolite) %>%
  group_modify(\(.x, .y) {
    a <- anova(lm(log2y ~ Lineage * Light, data = .x))
    tibble(p_interaction = a["Lineage:Light", "Pr(>F)"])
  }) %>%
  ungroup()

supp21 <- log2fc %>%
  left_join(interaction_p, by = c("Platform", "Metabolite")) %>%
  group_by(Platform) %>%
  mutate(FDR = p.adjust(p_interaction, method = "BH")) %>%
  ungroup() %>%
  mutate(
    MS = Platform,
    sig_class = case_when(
      FDR < 0.05 & log2FC_Wa > log2FC_Sp ~ "Wa > Sp",
      FDR < 0.05 & log2FC_Sp > log2FC_Wa ~ "Sp > Wa",
      TRUE ~ "Not significant"
    )
  ) %>%
  select(
    Platform, MS, Ion, Metabolite,
    log2LL_mean_Sp, log2HL_mean_Sp, log2FC_Sp,
    log2LL_mean_Wa, log2HL_mean_Wa, log2FC_Wa,
    p_interaction, FDR, sig_class
  ) %>%
  arrange(Platform, Metabolite)

write_csv(supp21, "output/SupplementaryData21_all_metabolite_log2FC_statistics.csv")

# Number of metabolites per platform and Lineage x Light interaction classes (FDR < 0.05)
# -> Supplementary Data 21: 1,730 metabolites
print(supp21 %>% count(Platform, sig_class) %>% pivot_wider(names_from = sig_class, values_from = n))
cat("Total metabolites:", nrow(supp21), "\n")
