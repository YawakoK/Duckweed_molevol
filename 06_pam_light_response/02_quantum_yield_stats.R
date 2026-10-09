# Step 02: quantum-yield partitioning [Y(II), Y(NPQ), Y(NO)] at two actinic light
#          (AL) intensities: summary, species contrasts and species x growth-light
#          interaction tests (Fig. 4d; Supplementary Data 15-17).
#
# Input : output/PAM_light_curves_processed.csv          (01)
# Output: output/SupplementaryData15_PAM_quantum_yield_summary.csv
#         output/SupplementaryData16_PAM_species_pairwise_contrasts_by_condition.csv
#         output/SupplementaryData17_PAM_species_growth_light_interaction_ANOVA.csv
#
# AL50   = measured PPFD 61.885 umol m-2 s-1 (light step 2)
# AL1000 = measured PPFD 1012.5 umol m-2 s-1 (light step 12)
# Models:
#   SD16: per yield, Y ~ species * growth_light * AL + (1 | sample_rep) (lme4/lmerTest),
#         emmeans pairs(species | growth_light * AL), Sidak adjustment
#         (one contrast per family, so the adjustment leaves P unchanged);
#         emmeans default d.f. for lmer models (Kenward-Roger, via pbkrtest).
#   SD17: per yield and AL, lm(Y ~ species * growth_light), sequential ANOVA,
#         species:growth_light term.
# Author: Natsu Katayama
# Run from this directory:  Rscript 02_quantum_yield_stats.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

OUT <- "output"

format_p <- function(p) {
  case_when(
    is.na(p) ~ "p=NA",
    p < 0.001 ~ "p<0.001",
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 3))
  )
}

growth_labels  <- c("PFD_50" = "LL", "PFD_1000" = "HL")
actinic_labels <- c("PPFD_50" = "AL50", "PPFD_1000" = "AL1000")

qy <- read.csv(file.path(OUT, "PAM_light_curves_processed.csv")) |>
  filter(abs(PPFD - 61.885) < 0.01 | abs(PPFD - 1012.500) < 0.01) |>
  mutate(
    PPFD_label = factor(if_else(abs(PPFD - 61.885) < 0.01, "PPFD_50", "PPFD_1000"),
                        levels = c("PPFD_50", "PPFD_1000")),
    pfd = factor(pfd, levels = c("PFD_50", "PFD_1000")),
    species = factor(recode(SpID, "Sp7498" = "Sp", "Wa8730" = "Wa"), levels = c("Sp", "Wa")),
    YII = Y.II., YNPQ = Y.NPQ., YNO = Y.NO.
  ) |>
  pivot_longer(c(YII, YNPQ, YNO), names_to = "Ytype", values_to = "Y") |>
  mutate(Ytype = factor(Ytype, levels = c("YII", "YNPQ", "YNO"),
                        labels = c("Y(II)", "Y(NPQ)", "Y(NO)")))

relabel <- function(d) {
  d |> mutate(
    growth_light  = factor(growth_labels[as.character(pfd)], levels = c("LL", "HL")),
    actinic_light = factor(actinic_labels[as.character(PPFD_label)], levels = c("AL50", "AL1000"))
  )
}

# ---- Supplementary Data 15: summary (n, mean, SD, SE) ----------------------
sd15 <- qy |>
  group_by(species, pfd, PPFD_label, Ytype) |>
  summarise(n = sum(!is.na(Y)), mean = mean(Y, na.rm = TRUE), sd = sd(Y, na.rm = TRUE),
            .groups = "drop") |>
  mutate(se = sd / sqrt(n)) |>
  arrange(species, pfd, PPFD_label, Ytype) |>
  relabel() |>
  select(species, growth_light, actinic_light, Ytype, n, mean, sd, se)
write.csv(sd15, file.path(OUT, "SupplementaryData15_PAM_quantum_yield_summary.csv"),
          row.names = FALSE)
# -> Methods "Statistical analyses": n per species x growth light x AL combination
cat("n per species x growth light x AL (plant samples):", unique(sd15$n), "\n")

# ---- Supplementary Data 16: Sp - Wa contrasts within growth light x AL -----
species_contrasts <- function(yt) {
  m <- lmer(Y ~ species * pfd * PPFD_label + (1 | sample_rep),
            data = filter(qy, Ytype == yt))
  emm <- emmeans(m, ~ species | pfd * PPFD_label)
  as.data.frame(pairs(emm, adjust = "sidak")) |> mutate(Ytype = yt)
}
sd16 <- bind_rows(lapply(levels(qy$Ytype), species_contrasts)) |>
  mutate(Ytype = factor(Ytype, levels = levels(qy$Ytype))) |>
  arrange(Ytype, PPFD_label, pfd) |>
  relabel() |>
  mutate(p_label = format_p(p.value)) |>
  select(Ytype, growth_light, actinic_light, contrast, estimate, SE, df, t.ratio,
         p.value, p_label)
write.csv(sd16, file.path(OUT, "SupplementaryData16_PAM_species_pairwise_contrasts_by_condition.csv"),
          row.names = FALSE)
# -> Results "Photosynthetic light-response curves and PSII energy partitioning":
#    Y(II) Wa > Sp at AL50 LL P = 0.005, AL50 HL P < 0.001, AL1000 LL and HL P < 0.001;
#    Y(NPQ) Sp > Wa at AL1000 in HL-grown plants P = 0.002
cat("\nSupplementary Data 16 (Sp - Wa, Sidak):\n")
print(sd16 |> select(Ytype, growth_light, actinic_light, estimate, df, t.ratio, p.value, p_label),
      digits = 4)

# ---- Supplementary Data 17: species x growth-light interaction per AL ------
interaction_term <- function(d) {
  a <- as.data.frame(anova(lm(Y ~ species * pfd, data = d)))
  a <- a["species:pfd", ]
  data.frame(sumsq = a$`Sum Sq`, meansq = a$`Mean Sq`, df = a$Df,
             F.value = a$`F value`, p.value = a$`Pr(>F)`)
}
sd17 <- qy |>
  group_by(PPFD_label, Ytype) |>
  group_modify(~ interaction_term(.x)) |>
  ungroup() |>
  mutate(actinic_light = factor(actinic_labels[as.character(PPFD_label)],
                                levels = c("AL50", "AL1000")),
         term = "species:growth_light",
         p_label = format_p(p.value)) |>
  arrange(actinic_light, Ytype) |>
  select(Ytype, actinic_light, term, sumsq, meansq, df, F.value, p.value, p_label)
write.csv(sd17, file.path(OUT, "SupplementaryData17_PAM_species_growth_light_interaction_ANOVA.csv"),
          row.names = FALSE)
# -> Results: no species x growth-light interaction for any yield at either AL (all P > 0.05)
cat("\nSupplementary Data 17 (species x growth light, lm per AL):\n")
print(sd17 |> select(Ytype, actinic_light, F.value, p.value), digits = 4)
cat("min interaction P =", signif(min(sd17$p.value), 3), "\n")
