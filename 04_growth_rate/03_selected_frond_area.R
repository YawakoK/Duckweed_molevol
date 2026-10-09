# Step 03: image-level frond area from the manually selected candidate mask.
#
# For every well image the frond area of the candidate mask (filter x preprocessing)
# chosen in the step-02 viewer is taken from the candidate-area table. The script checks
# that this equals the `selected_area` recorded by the viewer export.
#
# Input : data/image_mask_selection.csv                 manual mask choice per image (viewer export)
#         data/frond_area_candidate_masks_by_image.csv  12 candidate areas per image (step 02)
# Output: output/frond_area_selected_by_image.csv       720 images (90 wells x Day0-Day7)
# Run from this directory:  Rscript 03_selected_frond_area.R
# Author: Natsu Katayama
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

dir.create("output", showWarnings = FALSE)

selection <- read_csv("data/image_mask_selection.csv", show_col_types = FALSE,
                      col_types = cols(date = col_character(), .default = col_guess()))
candidates <- read_csv("data/frond_area_candidate_masks_by_image.csv", show_col_types = FALSE,
                       col_types = cols(date = col_character(), .default = col_guess()))

area <- selection |>
  select(id, selected_filter, selected_preprocess, selection_scope, selected_area) |>
  inner_join(candidates,
             by = c("id", "selected_filter" = "filter", "selected_preprocess" = "preprocess"))

stopifnot(nrow(area) == nrow(selection),
          max(abs(area$area_mm2 - area$selected_area)) < 1e-9)

area <- area |>
  transmute(id, batch, date, day, sample, light, wellID,
            sampleID2 = paste0(batch, "_", wellID), PhotoID, elapsed_time,
            selected_filter, selected_preprocess, selection_scope,
            frond_area = area_mm2)

write_csv(area, "output/frond_area_selected_by_image.csv")

cat("Images:", nrow(area), " wells:", n_distinct(area$sampleID2), "\n")
print(as.data.frame(count(area, selected_preprocess, selected_filter)))
