# Step 01: read the raw PAM light-curve exports (one file per
#          plant sample), attach the measured actinic PPFD of each light step,
#          compute ETR and write one long table for S. polyrhiza and W. australiana.
#
# Input : data/raw_PAM_light_curves/*.CSV      raw instrument exports (";"-separated)
#         data/raw_PAM_light_curves/PARlist_240702.csv
#                                              measured PPFD of each actinic light step
# Output: output/PAM_light_curves_processed.csv
#
# Growth light: PFD_50 = LL (growth PPFD 50), PFD_1000 = HL (growth PPFD 1000).
# ETR = Y(II) * PAR * 0.84 * 0.5, where PAR is the actinic PAR value recorded by
# the instrument for each step (column "PAR" of the raw export).
# Author: Natsu Katayama
# Run from this directory:  Rscript 01_process_light_curves.R
suppressPackageStartupMessages({
  library(dplyr)
})

RAW <- file.path("data", "raw_PAM_light_curves")
OUT <- "output"
dir.create(OUT, showWarnings = FALSE)

# Measured PPFD of the 17 actinic light steps of the light-curve program.
par_list <- read.csv(file.path(RAW, "PARlist_240702.csv"), fileEncoding = "UTF-8-BOM")

# One row per light-curve measurement (one plant sample each).
# Sp/Wa LL replicate 4 (measured 2024-07-05) is not used: the growth chamber
# light on that day was suspected to be brighter than the LL setting.
samples <- tribble(
  ~file,                                   ~sample,      ~rep,   ~SpID,    ~pfd,
  "Sp_50_rep1_240702_105813.CSV",          "Sp_pfd50",   "rep1", "Sp7498", "PFD_50",
  "Sp_50_rep2_240703_110054.CSV",          "Sp_pfd50",   "rep2", "Sp7498", "PFD_50",
  "Sp_50_rep3_240704_115929.CSV",          "Sp_pfd50",   "rep3", "Sp7498", "PFD_50",
  "Wa_50_rep1_20240702_122530.CSV",        "Wa_pfd50",   "rep1", "Wa8730", "PFD_50",
  "Wa_50_rep2_240703_122733.CSV",          "Wa_pfd50",   "rep2", "Wa8730", "PFD_50",
  "Wa_50_rep3_240704_131834.CSV",          "Wa_pfd50",   "rep3", "Wa8730", "PFD_50",
  "Sp_1000_DAY4_rep1_240712_100437.CSV",   "Sp_pfd1000", "rep1", "Sp7498", "PFD_1000",
  "Sp_1000_DAY4_rep2_240712_154305.CSV",   "Sp_pfd1000", "rep2", "Sp7498", "PFD_1000",
  "Sp_1000_LC_rep3_240725_105556.CSV",     "Sp_pfd1000", "rep3", "Sp7498", "PFD_1000",
  "Wa_1000_DAY4_rep1_240712_111242.CSV",   "Wa_pfd1000", "rep1", "Wa8730", "PFD_1000",
  "Wa_1000_DAY4_rep2_240712_132812.CSV",   "Wa_pfd1000", "rep2", "Wa8730", "PFD_1000",
  "Wa_1000_LC_rep3_240725_120403.CSV",     "Wa_pfd1000", "rep3", "Wa8730", "PFD_1000"
)

read_lc <- function(i) {
  s <- samples[i, ]
  read.delim(file.path(RAW, s$file), sep = ";", header = TRUE, na.strings = c(" ")) |>
    filter(!is.na(NPQ)) |>                       # keep the 17 light steps only
    mutate(sample = s$sample, rep = s$rep, SpID = s$SpID, pfd = s$pfd) |>
    cbind(par_list)
}

lc <- bind_rows(lapply(seq_len(nrow(samples)), read_lc)) |>
  mutate(
    ETR2 = Y.II. * PAR * 0.84 * 0.5,
    sample_rep = paste(sample, rep, sep = "_"),
    across(c(Fo., Fm.), as.numeric),
    SpID = factor(SpID, levels = c("Sp7498", "Wa8730")),
    pfd = factor(pfd, levels = c("PFD_50", "PFD_1000"))
  ) |>
  select(sample, rep, SpID, pfd, sample_rep, Setting, PPFD,
         F, Fo., Fm., X.Fo., Y.II., Y.NPQ., Y.NO., NPQ, qN, qP, qL, ETR, ETR2)

write.csv(lc, file.path(OUT, "PAM_light_curves_processed.csv"), row.names = FALSE)

cat("Light-curve measurements per species x growth light:\n")
print(lc |> distinct(SpID, pfd, sample_rep) |> count(SpID, pfd))
cat("Rows written:", nrow(lc), "\n")
