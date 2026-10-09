# Step 01: Read the per-sample metabolite abundance matrices of the three MS
#          platforms and reshape them into one long table.
#
# Input : data/CEMS_abundance.csv  (CE-TOF MS, 259 annotated features)
#         data/GCMS_abundance.csv  (GC-MS/MS, 167 annotated features)
#         data/LCMS_abundance.csv  (LC-MS/MS, BGI preprocessed table, 1,304 annotated features)
# Output: output/metabolite_abundance_long.csv
#         (Platform, Ion, Metabolite, sample, Lineage, Light, rep, value)
# Run from this directory: Rscript 01_read_abundance.R
# Author: Natsu Katayama
#
# Filtering applied (already reflected in the shipped matrices; re-applied here
# so that the rule is explicit): features without an identification
# ("Unknown" in CE-TOF MS / GC-MS/MS, "unidentified" in LC-MS/MS) are removed,
# and missing values of CE-TOF MS / GC-MS/MS are set to 0 (not detected).
# Samples: <Sp|Wa>_<HL|LL>_rep<1-5>; five biological replicates per group.
# LC-MS/MS pooled QC injections (columns "*_QC_*") are not used.

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

dir.create("output", showWarnings = FALSE)

sample_regex <- "^(Sp|Wa)_(HL|LL)_rep[0-9]+$"

read_platform <- function(file, platform, name_col, unknown_label, na_to_zero) {
  raw <- read.csv(file, check.names = FALSE)
  names(raw)[1] <- "Metabolite"
  raw <- raw %>% filter(.data[[name_col]] != unknown_label)
  if (!"Ion" %in% names(raw)) raw$Ion <- NA_character_
  samples <- grep(sample_regex, names(raw), value = TRUE)
  if (na_to_zero) raw[samples][is.na(raw[samples])] <- 0
  raw %>%
    select(Metabolite, Ion, all_of(samples)) %>%
    pivot_longer(all_of(samples), names_to = "sample", values_to = "value") %>%
    mutate(
      Platform = platform,
      Lineage = str_extract(sample, "^(Sp|Wa)"),
      Light = str_match(sample, "_(HL|LL)_")[, 2],
      rep = str_extract(sample, "rep[0-9]+")
    ) %>%
    select(Platform, Ion, Metabolite, sample, Lineage, Light, rep, value)
}

abundance_long <- bind_rows(
  read_platform("data/CEMS_abundance.csv", "CE-TOF MS", "Metabolite", "Unknown", TRUE),
  read_platform("data/GCMS_abundance.csv", "GC-MS/MS", "Metabolite", "Unknown", TRUE),
  read_platform("data/LCMS_abundance.csv", "LC-MS/MS", "Name", "unidentified", FALSE)
)

stopifnot(!any(is.na(abundance_long$value)))

n_features <- abundance_long %>% distinct(Platform, Metabolite) %>% count(Platform)
print(n_features)
# -> Supplementary Data 21: 259 (CE-TOF MS) + 167 (GC-MS/MS) + 1,304 (LC-MS/MS) = 1,730 metabolites
cat("Total metabolites:", sum(n_features$n), "\n")
print(abundance_long %>% distinct(Platform, sample, Lineage, Light) %>% count(Platform, Lineage, Light))

write_csv(abundance_long, "output/metabolite_abundance_long.csv")
