# Step 02: Fv/Fm statistics (Supplementary Data 11-14).
#
# Day 0, Day 1 and Day 4 only (6 lineages x 3 days x 6 = 108 measurements):
#   lm(Fv_Fm ~ sample * time)  ->  car::Anova(type = 2)            (Supplementary Data 13)
#   emmeans(~ time | sample), pairs(adjust = "tukey")                 (Supplementary Data 12)
#   emmeans(~ sample | time), multcomp::cld(adjust = "sidak")         (Supplementary Data 14)
# The compact letters of Supplementary Data 14 / Supplementary Fig. 3 were computed with
# Sidak-adjusted pairwise comparisons (and Sidak-adjusted 95% CIs); with Tukey adjustment the
# letters are identical (checked at the end of this script).
#
# Input : output/SupplementaryData10_FvFm_extracted_data.csv   (step 01)
# Output: output/SupplementaryData11_FvFm_summary_by_species_and_time.csv
#         output/SupplementaryData12_FvFm_time_pairwise_comparisons_within_species.csv
#         output/SupplementaryData13_FvFm_species_time_ANOVA.csv
#         output/SupplementaryData14_FvFm_species_group_letters_within_time.csv
# Run from this directory:  Rscript 02_fvfm_stats.R
# Author: Natsu Katayama

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(car)
  library(emmeans)
  library(multcomp)
})

OUT <- "output"
analysis_samples <- c("Sp7498", "Lp9387", "Lgp8L", "LaM10E", "Wh9525", "Wa")
analysis_times <- c("DAY0", "DAY1", "DAY4")

fvfm <- read_csv(file.path(OUT, "SupplementaryData10_FvFm_extracted_data.csv"), show_col_types = FALSE) |>
  filter(time %in% analysis_times) |>
  mutate(
    sample = factor(sample, levels = analysis_samples),
    time = factor(time, levels = analysis_times),
    antho = factor(antho, levels = c("antho", "non_antho"))
  )

# ---- Supplementary Data 11: mean, SD, SE ----
summary_df <- fvfm |>
  group_by(time, sample, antho) |>
  summarise(
    mean_FvFm = mean(Fv_Fm, na.rm = TRUE),
    sd_FvFm = sd(Fv_Fm, na.rm = TRUE),
    n = sum(!is.na(Fv_Fm)),
    se_FvFm = sd_FvFm / sqrt(n),
    .groups = "drop"
  )
write.csv(summary_df, file.path(OUT, "SupplementaryData11_FvFm_summary_by_species_and_time.csv"), row.names = FALSE)

# ---- Supplementary Data 13: linear model, Type II ANOVA ----
m_fvfm <- lm(Fv_Fm ~ sample * time, data = fvfm)
anova_res <- as.data.frame(car::Anova(m_fvfm, type = 2))
anova_res <- cbind(Effect = rownames(anova_res), anova_res)
write.csv(anova_res, file.path(OUT, "SupplementaryData13_FvFm_species_time_ANOVA.csv"), row.names = FALSE, na = "")

# ---- Supplementary Data 12: Day contrasts within each lineage (Tukey) ----
time_pairs <- pairs(emmeans(m_fvfm, ~ time | sample), adjust = "tukey") |>
  as.data.frame() |>
  dplyr::select(contrast, estimate, SE, df, t.ratio, p.value, sample)
write.csv(time_pairs, file.path(OUT, "SupplementaryData12_FvFm_time_pairwise_comparisons_within_species.csv"), row.names = FALSE)

# ---- Supplementary Data 14: compact letters among lineages within each day (Sidak) ----
emm_sample_by_time <- emmeans(m_fvfm, ~ sample | time)
cld_res <- multcomp::cld(emm_sample_by_time, adjust = "sidak", Letters = letters) |>
  as.data.frame()
cld_res$.group <- gsub(" ", "", cld_res$.group)
cld_res <- cld_res |> dplyr::select(time, SE, df, sample, emmean, lower.CL, upper.CL, .group)
write.csv(cld_res, file.path(OUT, "SupplementaryData14_FvFm_species_group_letters_within_time.csv"), row.names = FALSE)

# ---- values cited in the main text ----
cat("\nType II ANOVA (Supplementary Data 13)\n")
print(anova_res, row.names = FALSE)
# -> Results "Lineage-specific photoinhibition and recovery under high light":
#    lineage, time and lineage x time all P < 0.001
#    (sample F5,90 = 48.76, P = 3.6e-24; time F2,90 = 569.0, P = 8.4e-52; interaction F10,90 = 12.99, P = 9.5e-14)

cat("\nDay 1 vs Day 4 within lineage (Tukey; Supplementary Data 12, Fig. 4b)\n")
print(time_pairs |> filter(contrast == "DAY1 - DAY4") |> dplyr::select(sample, estimate, t.ratio, p.value), row.names = FALSE)
# -> Results: significant increase Day 1 -> Day 4 in S. polyrhiza (P = 1.7e-08), L. gibba (P = 0.012),
#    W. australiana (P = 0.0068); not in L. punctata (P = 0.24), L. aequinoctialis (P = 0.77),
#    Wolffiella hyalina (P = 0.93)

cat("\nHighest / lowest mean Fv/Fm per day (Supplementary Data 11)\n")
print(summary_df |> group_by(time) |>
        summarise(highest = sample[which.max(mean_FvFm)], second = sample[order(-mean_FvFm)[2]],
                  lowest = sample[which.min(mean_FvFm)], second_lowest = sample[order(mean_FvFm)[2]]))
# -> Results: Day 0 highest W. australiana, lowest Wolffiella hyalina;
#    Day 1 highest W. australiana and Wolffiella hyalina, lowest S. polyrhiza and L. gibba;
#    Day 4 highest W. australiana, then Wolffiella hyalina, lowest L. gibba

# Adjustment check: letters with Tukey-adjusted instead of Sidak-adjusted pairwise tests
cld_tukey <- suppressMessages(multcomp::cld(emm_sample_by_time, adjust = "tukey", Letters = letters)) |>
  as.data.frame()
cld_tukey$.group <- gsub(" ", "", cld_tukey$.group)
chk <- merge(cld_res[, c("time", "sample", ".group")], cld_tukey[, c("time", "sample", ".group")],
             by = c("time", "sample"))
cat("\nSidak and Tukey compact letters identical:", identical(chk$.group.x, chk$.group.y), "\n")
