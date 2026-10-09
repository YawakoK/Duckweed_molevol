# Step 03: species matching, coordinate cleaning and 0.5-degree spatial thinning.
#
#   1. drop records without species-level identification or coordinates
#   2. binomial = first two words of scientificName; keep the 36 accepted species
#      of data/species_list_anthocyanin_status.tsv and attach anthocyanin status
#   3. drop out-of-range coordinates, uncertainty > 10 km (missing kept), (0, 0)
#   4. drop duplicated species x coordinate records (first kept)
#   5. keep one random record per species per 0.5 x 0.5 degree cell
#      (cell = round(coord * 2) / 2; set.seed(123))
# The thinned set is checked against Supplementary Data 2 (gbifID, same order).
#
# Input : data/gbif_occurrences_compact.csv.xz (step 02)
#         data/species_list_anthocyanin_status.tsv
#         data/SupplementaryData02_GBIF_thinned_occurrence_records.csv (for the check)
# Output: output/filtering_counts.csv
#         output/thinned_records.csv
# Run from this directory:  Rscript 03_filter_and_thin.R
suppressPackageStartupMessages(library(tidyverse))

raw <- read_csv("data/gbif_occurrences_compact.csv.xz", na = "NA",
                col_types = cols(source_file = "c", gbifID = "c", species = "c",
                                 scientificName = "c", decimalLatitude = "d",
                                 decimalLongitude = "d",
                                 coordinateUncertaintyInMeters = "d"))

# Species list (trailing spaces in two names are removed so that they join)
sp_list <- read.delim("data/species_list_anthocyanin_status.tsv", header = TRUE, sep = "\t") |>
  mutate(across(where(is.character), str_squish))

gbif <- raw |>
  filter(species != "") |>
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) |>
  mutate(spName = str_extract(scientificName, "^[^ ]+ [^ ]+"),
         author = sub("^[^ ]+ [^ ]+ ", "", scientificName))

matched <- gbif |>
  filter(spName %in% sp_list$sp) |>
  left_join(sp_list |> select(genus_cor, sp, antho), by = c("spName" = "sp"))

cat("Accepted species in the list:", nrow(sp_list), "\n")   # -> Methods: 36 species
cat("Species without GBIF records:",
    paste(setdiff(sp_list$sp, matched$spName), collapse = ", "), "\n")

clean0 <- matched |>
  filter(decimalLatitude >= -90, decimalLatitude <= 90,
         decimalLongitude >= -180, decimalLongitude <= 180) |>
  mutate(antho = factor(antho, levels = c("Antho", "NonAntho")))
clean1 <- clean0 |>
  filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 10000)
clean2 <- clean1 |>
  filter(!(decimalLatitude == 0 & decimalLongitude == 0))
clean3 <- clean2 |>
  distinct(spName, decimalLatitude, decimalLongitude, .keep_all = TRUE)

set.seed(123)
thinned <- clean3 |>
  mutate(lat_bin = round(decimalLatitude * 2) / 2,
         lon_bin = round(decimalLongitude * 2) / 2) |>
  group_by(spName, lat_bin, lon_bin) |>
  slice_sample(n = 1) |>
  ungroup()

counts <- tibble(
  step = c("raw GBIF records",
           "species-level identification and coordinates present",
           "assigned to an accepted species",
           "valid coordinate range",
           "coordinate uncertainty <= 10 km or missing",
           "not at (0, 0)",
           "species x coordinate duplicates removed",
           "0.5-degree spatial thinning"),
  n_records = c(nrow(raw), nrow(gbif), nrow(matched), nrow(clean0), nrow(clean1),
                nrow(clean2), nrow(clean3), nrow(thinned)),
  n_species = c(NA, NA, n_distinct(matched$spName), n_distinct(clean0$spName),
                n_distinct(clean1$spName), n_distinct(clean2$spName),
                n_distinct(clean3$spName), n_distinct(thinned$spName)))
print(counts, width = Inf)
# -> Results / Methods: 18,227 occurrence records from 29 species (last row)
write_csv(counts, "output/filtering_counts.csv")

thinned_out <- thinned |>
  mutate(abs_latitude = abs(decimalLatitude)) |>
  select(source_file, gbifID, spName, genus_cor, antho, scientificName,
         decimalLatitude, decimalLongitude, abs_latitude, lat_bin, lon_bin,
         coordinateUncertaintyInMeters)
write_csv(thinned_out, "output/thinned_records.csv", na = "NA")

# Check against the submitted Supplementary Data 2
sd2 <- read_csv("data/SupplementaryData02_GBIF_thinned_occurrence_records.csv",
                col_types = cols(gbifID = "c", .default = "?"))
cat("Identical to Supplementary Data 2 (gbifID, row order):",
    identical(thinned_out$gbifID, sd2$gbifID), "\n")
cat("Identical coordinates and species:",
    isTRUE(all.equal(thinned_out$decimalLatitude, sd2$decimalLatitude)) &&
    isTRUE(all.equal(thinned_out$decimalLongitude, sd2$decimalLongitude)) &&
    identical(thinned_out$spName, sd2$spName) &&
    identical(as.character(thinned_out$antho), sd2$antho), "\n")
