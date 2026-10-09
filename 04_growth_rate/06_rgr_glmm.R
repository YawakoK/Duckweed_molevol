# Step 06: Gamma GLMMs for Day 0 - Day 4 RGR (Fig. 3; Supplementary Data 7-9).
#
# Models (lme4::glmer, Gamma(link = "log"), nAGQ = 0, bobyqa):
#   species model     : RGR ~ species * light + (1 | batch)                     (Fig. 3a)
#   anthocyanin model : RGR ~ antho * light + (1 | species) + (1 | batch)       (Fig. 3b)
# Type II Wald chi-square tests: car::Anova(type = 2).
# emmeans on the response scale: HL/LL ratio within each species (Fig. 3a labels) and
# non-antho/antho ratio within each light condition (Fig. 3b), unadjusted.
#
# Input : output/SupplementaryData06_RGR_by_well_Day0_Day4.csv (04)
# Output: output/SupplementaryData07_RGR_anthocyanin_light_TypeII_Wald.csv
#         output/SupplementaryData08_RGR_EMMs_by_anthocyanin_status_and_light.csv
#         output/SupplementaryData09_RGR_pairwise_anthocyanin_contrasts_by_light.csv
#         output/Fig3a_species_light_TypeII_Wald.csv
#         output/Fig3a_HL_LL_ratio_by_species.csv
# Run from this directory:  Rscript 06_rgr_glmm.R
# Author: Natsu Katayama
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(lme4)
  library(car)
  library(emmeans)
})

fig3_samples <- c("Sp", "Lp", "Lgp8L", "M10E", "Wh", "Wa")
light_levels <- c("PFD100", "PFD1000")   # PFD100 = LL, PFD1000 = HL

rgr <- read_csv("output/SupplementaryData06_RGR_by_well_Day0_Day4.csv", show_col_types = FALSE) |>
  mutate(
    sample = factor(sample, levels = fig3_samples),
    light = factor(light, levels = light_levels),
    antho = factor(antho, levels = c("antho", "non-antho"))
  )

fit_gamma_glmm <- function(formula, data) {
  glmer(formula, data = data, family = Gamma(link = "log"), nAGQ = 0,
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}

## species model (Fig. 3a) --------------------------------------------------------------------
model_sample <- fit_gamma_glmm(rgr_per_day ~ sample * light + (1 | batch), rgr)
anova_sample <- as.data.frame(Anova(model_sample, type = 2)) |> rownames_to_column("Effect")
write.csv(anova_sample, "output/Fig3a_species_light_TypeII_Wald.csv", row.names = FALSE)

light_contrasts <- as.data.frame(
  emmeans(model_sample, pairwise ~ light | sample, type = "response")$contrasts
) |>
  mutate(hl_vs_ll_ratio = 1 / ratio)
write.csv(light_contrasts, "output/Fig3a_HL_LL_ratio_by_species.csv", row.names = FALSE)

## anthocyanin model (Fig. 3b) ----------------------------------------------------------------
model_antho <- fit_gamma_glmm(rgr_per_day ~ antho * light + (1 | sample) + (1 | batch), rgr)
anova_antho <- as.data.frame(Anova(model_antho, type = 2)) |> rownames_to_column("Effect")
write.csv(anova_antho, "output/SupplementaryData07_RGR_anthocyanin_light_TypeII_Wald.csv", row.names = FALSE)

emm_antho <- emmeans(model_antho, ~ antho | light, type = "response")
emm_antho_df <- as.data.frame(emm_antho) |>
  transmute(light, antho, estimated_RGR = response, SE, df,
            asymp_LCL = asymp.LCL, asymp_UCL = asymp.UCL)
write.csv(emm_antho_df, "output/SupplementaryData08_RGR_EMMs_by_anthocyanin_status_and_light.csv",
          row.names = FALSE)

pairwise_antho <- as.data.frame(contrast(emm_antho, method = "revpairwise", adjust = "none")) |>
  select(light, contrast, ratio, SE, df, null, z.ratio, p.value)
write.csv(pairwise_antho, "output/SupplementaryData09_RGR_pairwise_anthocyanin_contrasts_by_light.csv",
          row.names = FALSE)

## values cited in the main text --------------------------------------------------------------
cat("\nSpecies x light (Fig. 3a):\n"); print(anova_sample)
cat("\nHL/LL RGR ratio by species:\n")
print(light_contrasts |> select(sample, hl_vs_ll_ratio, p.value))
# -> Results "Growth rate comparisons": HL/LL ratio 1.29 (L. punctata, Lp) to 1.85 (W. hyalina, Wh)
cat("\nAnthocyanin status x light (Supplementary Data 7):\n"); print(anova_antho)
# -> Results: interaction P = 0.048 (row "antho:light")
cat("\nEstimated marginal means (Supplementary Data 8):\n"); print(emm_antho_df)
cat("\nnon-antho / antho ratio by light (Supplementary Data 9):\n"); print(pairwise_antho)
# -> Results: ratio 1.37 under LL (PFD100) and 1.57 under HL (PFD1000), both P <= 0.001
