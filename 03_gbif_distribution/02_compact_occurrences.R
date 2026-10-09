# Step 02: combine the five per-genus GBIF extracts into one compact table.
#
# Keeps every raw record (753,829) in the original row order (Spirodela, Landoltia,
# Lemna, Wolffiella, Wolffia), but only the columns needed for filtering and
# thinning, so that step 03 can be run from data/ without the 380 MB GBIF dump.
# Row order matters: step 03 keeps the first of duplicated coordinates and
# samples within grid cells with a fixed seed.
#
# Input : gbif_raw/<Genus>_extracted.txt  (step 01; not shipped)
# Output: data/gbif_occurrences_compact.csv.xz
# Run from this directory:  Rscript 02_compact_occurrences.R
suppressPackageStartupMessages(library(tidyverse))

RAW <- "gbif_raw"
genera <- c("Spirodela", "Landoltia", "Lemna", "Wolffiella", "Wolffia")

raw <- bind_rows(lapply(genera, function(g) {
  f <- paste0(g, "_extracted.txt")
  d <- read.delim(file.path(RAW, f), header = TRUE, sep = "\t", quote = "",
                  colClasses = c(gbifID = "character"))
  d$source_file <- f
  d
}))
cat("raw GBIF records:", nrow(raw), "\n")   # -> 753,829

compact <- raw |>
  select(source_file, gbifID, species, scientificName,
         decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters)

write_csv(compact, "data/gbif_occurrences_compact.csv.xz", na = "NA")
