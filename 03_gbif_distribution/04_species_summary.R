# Step 04: species-level latitude summary (Supplementary Data 1) and the
# median absolute latitudes by anthocyanin status cited in the main text.
#
# Input : data/SupplementaryData02_GBIF_thinned_occurrence_records.csv
#         (= Supplementary Data 2; reproduced by step 03)
# Output: output/SupplementaryData01_GBIF_species_and_latitude_summary.csv
#         output/abs_latitude_by_anthocyanin_status.csv
# Run from this directory:  Rscript 04_species_summary.R
suppressPackageStartupMessages(library(tidyverse))

records <- read.csv("data/SupplementaryData02_GBIF_thinned_occurrence_records.csv") |>
  mutate(antho = factor(antho, levels = c("Antho", "NonAntho")))

cat("Thinned records:", nrow(records), "  species:", n_distinct(records$spName), "\n")
# -> Results "...18,227 occurrence records from 29 species"

species_summary <- records |>
  group_by(spName, genus = genus_cor, antho) |>
  summarise(n_records = n(),
            median_abs_latitude = median(abs_latitude),
            mean_abs_latitude = mean(abs_latitude),
            min_latitude = min(decimalLatitude),
            max_latitude = max(decimalLatitude),
            .groups = "drop") |>
  arrange(antho, spName)
write.csv(species_summary,
          "output/SupplementaryData01_GBIF_species_and_latitude_summary.csv",
          row.names = FALSE)

group_summary <- bind_rows(
  records |> mutate(scope = "All genera"),
  records |> filter(genus_cor == "Lemna") |> mutate(scope = "Lemna")) |>
  group_by(scope, antho) |>
  summarise(n_records = n(),
            n_species = n_distinct(spName),
            median_abs_latitude = median(abs_latitude),
            mean_abs_latitude = mean(abs_latitude),
            .groups = "drop")
print(group_summary, width = Inf)
# -> Results: median absolute latitude 48.3 (accumulating) vs 37.5 (non-accumulating)
#    (rows "All genera"; Fig. 2b).  Rows "Lemna": Supplementary Fig. 1b.
write.csv(group_summary, "output/abs_latitude_by_anthocyanin_status.csv", row.names = FALSE)
