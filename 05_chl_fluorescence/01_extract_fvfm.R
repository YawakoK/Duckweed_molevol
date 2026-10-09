# Step 01: extract the Fv/Fm measurements used for Fig. 4a-b and Supplementary Fig. 3
#          from the raw PAM export tables (Supplementary Data 10).
#
# Extraction rules
#   - growth / exposure light pfd == 1000 (HL), time points DAY0-DAY4
#   - S. polyrhiza (Sp7498), L. punctata (Lp9387), L. aequinoctialis (LaM10E),
#     W. hyalina (Wh9525) and W. australiana (Wa) from the nine dated tables
#     dat240408 ... dat250609
#   - L. gibba (Lgp8L) from dat1-dat5 only (time recorded as "1d3h", "2d", ...:
#     converted to elapsed days; the 3 h point (0.125 d) is excluded, the rest
#     rounded to DAY0-DAY4)
#   - the last six measurements (file order) are kept for every lineage x time
#   -> 6 lineages x 5 days x 6 = 180 measurements
#
# Input : data/dat240408.txt, dat240520.txt, dat240527.txt, dat240610.txt,
#         dat240617.txt, dat240624.txt, dat240722.txt, dat250420.txt,
#         dat250609.txt, dat1.txt - dat5.txt   (raw PAM exports, tab-separated)
# Output: output/SupplementaryData10_FvFm_extracted_data.csv
# Run from this directory:  Rscript 01_extract_fvfm.R
# Author: Natsu Katayama

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

DATA <- "data"
OUT <- "output"
dir.create(OUT, showWarnings = FALSE)

extract_times <- paste0("DAY", 0:4)
current_samples <- c("Sp7498", "Lp9387", "LaM10E", "Wh9525", "Wa")
analysis_samples <- c("Sp7498", "Lp9387", "Lgp8L", "LaM10E", "Wh9525", "Wa")
antho_samples <- c("Sp7498", "Lp9387", "Lgp8L")   # anthocyanin-accumulating lineages

# "DAY4" -> 4, "1d3h" -> 1.125, "7d" -> 7
parse_elapsed_days <- function(x) {
  x <- as.character(x)
  elapsed <- rep(NA_real_, length(x))
  day_idx <- grepl("^DAY[0-9]+$", x)
  elapsed[day_idx] <- as.numeric(sub("^DAY", "", x[day_idx]))
  dh_idx <- grepl("^[0-9]+d[0-9]+h$", x)
  elapsed[dh_idx] <- as.numeric(sub("d.*$", "", x[dh_idx])) +
    as.numeric(sub("^.*d([0-9]+)h$", "\\1", x[dh_idx])) / 24
  d_idx <- grepl("^[0-9]+d$", x)
  elapsed[d_idx] <- as.numeric(sub("d$", "", x[d_idx]))
  elapsed
}

read_fvfm <- function(file_name) {
  read.delim(
    file.path(DATA, file_name),
    sep = "\t", header = TRUE,
    na.strings = c("", " ", "na", "NA"),
    check.names = FALSE, colClasses = "character"
  ) |>
    mutate(source_file = file_name)
}

# The file order matters: slice_tail() keeps the last six rows per lineage x time.
current_file_names <- c(
  "dat240408.txt", "dat240527.txt", "dat240520.txt", "dat240610.txt",
  "dat240617.txt", "dat240624.txt", "dat240722.txt", "dat250420.txt",
  "dat250609.txt"
)

# as.numeric() warns "NAs introduced by coercion" for one Fv_Fm cell recorded as "-"
# in dat240520.txt (not a selected measurement).
current_5sp <- bind_rows(lapply(current_file_names, read_fvfm)) |>
  mutate(
    sample = recode(sample, "SP7498" = "Sp7498"),
    across(c(F0, Fm, Fv_Fm, pfd), as.numeric)
  ) |>
  filter(pfd == 1000, time %in% extract_times, sample %in% current_samples) |>
  group_by(sample, time) |>
  slice_tail(n = 6) |>
  ungroup() |>
  select(sample, F0, Fm, Fv_Fm, pfd, time, ld, source_file)

lgp8l <- bind_rows(lapply(paste0("dat", 1:5, ".txt"), read_fvfm)) |>
  mutate(
    sample = recode(sample, "2_Lgp8L" = "Lgp8L"),
    across(c(F0, Fm, Fv_Fm, pfd), as.numeric),
    elapsed_days = parse_elapsed_days(time),
    time = paste0("DAY", round(elapsed_days))
  ) |>
  filter(
    sample == "Lgp8L", pfd == 1000,
    elapsed_days != 0.125,            # drop the 3 h time point
    time %in% extract_times,
    is.finite(Fv_Fm)
  ) |>
  group_by(sample, time) |>
  slice_tail(n = 6) |>
  ungroup() |>
  select(sample, F0, Fm, Fv_Fm, pfd, time, ld, source_file)

fvfm <- bind_rows(current_5sp, lgp8l) |>
  mutate(
    sample = factor(sample, levels = analysis_samples),
    time = factor(time, levels = extract_times),
    antho = if_else(as.character(sample) %in% antho_samples, "antho", "non_antho")
  ) |>
  arrange(sample, time, source_file) |>
  select(sample, F0, Fm, Fv_Fm, pfd, time, ld, antho)

write_csv(fvfm, file.path(OUT, "SupplementaryData10_FvFm_extracted_data.csv"))

cat("Supplementary Data 10:", nrow(fvfm), "measurements\n")   # -> Methods "Statistical analysis": n = 180 measurements
print(table(fvfm$sample, fvfm$time))                          # 6 per lineage x day
