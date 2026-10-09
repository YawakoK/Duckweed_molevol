# Step 04: Day 0 - Day 4 relative growth rate (RGR) per well (Supplementary Data 6).
#
#   RGR = [ln(A_Day4) - ln(A_Day0)] / (t_Day4 - t_Day0),  t in days from the first photo of the batch
# Wells with positive RGR are kept; for each species x light combination the six wells
# with the smallest Day 0 frond area are used (to limit early crowding), giving n = 72 wells.
# Anthocyanin status: Sp, Lp, Lgp8L = antho (accumulating); M10E, Wh, Wa = non-antho.
#
# Input : output/frond_area_selected_by_image.csv (03)
# Output: output/SupplementaryData06_RGR_by_well_Day0_Day4.csv
# Run from this directory:  Rscript 04_rgr_day0_day4.R
# Author: Natsu Katayama
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

fig3_samples <- c("Sp", "Lp", "Lgp8L", "M10E", "Wh", "Wa")
light_levels <- c("PFD100", "PFD1000")
antho_samples <- c("Sp", "Lp", "Lgp8L")

area <- read_csv("output/frond_area_selected_by_image.csv", show_col_types = FALSE) |>
  mutate(
    sample = factor(sample, levels = fig3_samples),
    light = factor(light, levels = light_levels),
    elapsed_days = elapsed_time / 24,
    antho = factor(if_else(as.character(sample) %in% antho_samples, "antho", "non-antho"),
                   levels = c("antho", "non-antho"))
  )

start_df <- area |>
  filter(day == "Day0", frond_area > 0) |>
  transmute(batch, wellID, sampleID2, sample, light, antho,
            start_day = day, start_elapsed_time = elapsed_time,
            start_elapsed_days = elapsed_days, start_frond_area = frond_area)

end_df <- area |>
  filter(day == "Day4", frond_area > 0) |>
  transmute(batch, wellID, sampleID2,
            end_day = day, end_elapsed_time = elapsed_time,
            end_elapsed_days = elapsed_days, end_frond_area = frond_area)

rgr_all <- start_df |>
  inner_join(end_df, by = c("batch", "wellID", "sampleID2")) |>
  mutate(
    window = "day0_day4",
    duration_days = end_elapsed_days - start_elapsed_days,
    rgr_per_day = (log(end_frond_area) - log(start_frond_area)) / duration_days,
    fold_change = end_frond_area / start_frond_area
  ) |>
  filter(is.finite(rgr_per_day), rgr_per_day > 0, duration_days > 0)

rgr <- rgr_all |>
  group_by(sample, light) |>
  arrange(start_frond_area, desc(rgr_per_day), .by_group = TRUE) |>
  slice_head(n = 6) |>
  ungroup()

write_csv(rgr, "output/SupplementaryData06_RGR_by_well_Day0_Day4.csv")

cat("Candidate wells:", nrow(rgr_all), " selected wells:", nrow(rgr), "\n")  # -> Methods: n = 72 wells
print(as.data.frame(count(rgr, sample, light)))
