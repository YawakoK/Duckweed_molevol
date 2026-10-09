# 04_growth_rate — Relative growth rate under low and high light

Produces: Fig. 3a, 3b; Supplementary Fig. 2; Supplementary Data 4, 5, 6, 7, 8, 9
(Supplementary Data 3, species and strain information, is included in `data/` for reference)
Author: Natsu Katayama

## Overview

Six species (Sp, Lp, Lgp8L, M10E, Wh, Wa) were grown in 6-well plates under LL (PFD100) and
HL (PFD1000) and photographed every ~24 h for seven days (five experimental batches).
Frond area was segmented from the photos with OpenCV (steps 01–02, manual choice of the
best candidate mask per image), and the R steps compute Day 0–Day 4 RGR, growth-curve fits and
the Gamma GLMMs.

| Step | Script | What it does | Output |
|---|---|---|---|
| 01 | `01_segment_wells_opencv.py` | ArUco plate cropping (DICT_4X4_50, 3000 × 2000 px), split into 2 × 3 wells, primary CIELAB mask (L 0–255, a 90–135, b 155–255) + thin-line/crescent removal; area = pixels × 0.0016 mm² | `processed/` images, `data/frond_area_main_mask_by_well_image.csv` |
| 02 | `02_candidate_masks_viewer.py` | 6 candidate masks × {original, CLAHE-brightened} image per well image (broader mask L 36–255, a 0–138, b 143–255; red/orange OR mask L 0–255, a 136–255, b 134–244; small-component, reflection and root-line removal); HTML viewer used for the manual mask selection | `data/frond_area_candidate_masks_by_image.csv`, `viewer/index.html` |
| — | (manual) | best-matching mask chosen per well / per image in the viewer, exported with "Export CSV" | `data/image_mask_selection.csv` |
| 03 | `03_selected_frond_area.R` | image-level frond area of the selected mask (checked against the viewer export) | `output/frond_area_selected_by_image.csv` |
| 04 | `04_rgr_day0_day4.R` | RGR = [ln A(Day4) − ln A(Day0)] / Δt (days); 6 wells per species × light with the smallest Day 0 area (n = 72) | `output/SupplementaryData06_RGR_by_well_Day0_Day4.csv` |
| 05 | `05_growth_curve_models.R` | exponential (lm on ln area) and logistic (nls) fits of Day 0–7 per well; AIC on log scale; exponential-fixed / logistic-fixed / AIC-selected strategies, each followed by the Gamma GLMMs | `output/SupplementaryData04_…csv`, `output/SupplementaryData05_…csv`, `output/SupplementaryFig2.pdf` |
| 06 | `06_rgr_glmm.R` | `glmer(RGR ~ species * light + (1\|batch))` and `glmer(RGR ~ antho * light + (1\|species) + (1\|batch))`, Gamma(log), nAGQ = 0; `car::Anova` type II; `emmeans` (response scale) | `output/SupplementaryData07_…csv`, `08_…csv`, `09_…csv`, `output/Fig3a_*.csv` |
| 07 | `07_fig3.R` | Fig. 3a, 3b | `output/Fig3.pdf` |

## Inputs (`data/`)

| File | Content | Source |
|---|---|---|
| `image_metadata/photoinfo_<date>.csv` | EXIF capture time and plate ID of every photo (5 batches; `<date>` = experiment start, DDMMYYYY) | exiftool export + plate annotation |
| `image_metadata/sample_list_<date>.csv` | plate/well → species, light, temperature | experiment records |
| `frond_area_main_mask_by_well_image.csv` | step-01 output: 1,008 well images (all species photographed, incl. species not used in Fig. 3), primary-mask area | step 01 (as run 2026-07-14; image paths made relative) |
| `frond_area_candidate_masks_by_image.csv` | step-02 output: 720 well images of the six species × 2 preprocessing × 6 masks = 8,640 candidate areas (mm²) | step 02 (as run 2026-07-14; extracted from the viewer data) |
| `image_mask_selection.csv` | manual selection per image: `selected_filter`, `selected_preprocess`, `selected_area`, `selection_scope` (`well` = well default, `image` = per-image override) | viewer export, unchanged |
| `SupplementaryData03_species_and_strain_information.csv` | = Supplementary Data 3 | — |

Raw photos (`photos/<date>/Day<N>/*.JPG`) and the intermediate cropped/binary images
(`processed/`) are deposited in the Zenodo dataset and are not part of this repository; steps 01–02
need them and OpenCV, all later steps run from `data/` as is.

Codes: `Sp` *Spirodela polyrhiza*, `Lp` *Landoltia punctata*, `Lgp8L` *Lemna gibba* (Lg),
`M10E` *Lemna aequinoctialis* (La), `Wh` *Wolffiella hyalina*, `Wa` *Wolffia australiana*;
`PFD100` = LL, `PFD1000` = HL; `antho` = anthocyanin-accumulating (Sp, Lp, Lgp8L),
`non-antho` = non-accumulating (M10E, Wh, Wa).

## How to run

```bash
# optional, needs photos/ from Zenodo and Python with opencv-python, numpy, pandas, pillow
python 01_segment_wells_opencv.py
python 02_candidate_masks_viewer.py      # then select masks in viewer/index.html and export
# from the shipped data/
Rscript 03_selected_frond_area.R         # each R step runs in < 5 s
Rscript 04_rgr_day0_day4.R               # Supplementary Data 6
Rscript 05_growth_curve_models.R         # Supplementary Data 4, 5; Supplementary Fig. 2
Rscript 06_rgr_glmm.R                    # Supplementary Data 7-9, Fig. 3 statistics
Rscript 07_fig3.R                        # Fig. 3
```

`set.seed(1)` is used only for the point jitter in Fig. 3. The submitted Fig. 3 was finished in
a vector editor (species abbreviations on the axes, panel letters, layout); data and statistical
annotations are those of `output/Fig3.pdf`. Supplementary Fig. 2 pages 1–3 are panels a–c.

## Values reported in the main text and where they are printed

| Main text | Value | Printed by |
|---|---|---|
| Results "Growth rate comparisons"; Methods | n = 72 wells (6 per species × light) | `04` |
| Results | HL/LL RGR ratio 1.29 (*L. punctata*) to 1.85 (*W. hyalina*) | `06` (`HL/LL RGR ratio by species`; `output/Fig3a_HL_LL_ratio_by_species.csv`, column `hl_vs_ll_ratio`) |
| Results | anthocyanin status × light interaction P = 0.048 | `06`; Supplementary Data 7 |
| Results | non-antho/antho RGR ratio 1.37 (LL), 1.57 (HL), both P ≤ 0.001 | `06`; Supplementary Data 9 (EMMs in Supplementary Data 8) |
| Fig. 3a | species × light P = 0.01 | `06` (`Species x light`); `output/Fig3a_species_light_TypeII_Wald.csv` |

## Verification

All scripts were re-run from a clean copy under R 4.6.0 (the manuscript analyses used R 4.4.1).
Supplementary Data 6 and 8 are byte-identical to the submitted files; Supplementary Data 4 and 5
agree to < 1e-14 (relative), Supplementary Data 7 and 9 to < 3e-8 (last digits of the Wald chi-square,
z and P values, attributable to the newer lme4 version). All main-text values are reproduced.

## Environment

R 4.6.0; lme4 2.0.1, car 3.1-5, emmeans 2.0.3, dplyr 1.2.1, readr 2.2.0, ggplot2 4.0.3,
patchwork 1.3.2. Full list: `../environment/R_sessionInfo_04_growth_rate.txt`.
Python 3 with opencv-python (cv2.aruco, ≥ 4.7), numpy, pandas, pillow for steps 01–02.
