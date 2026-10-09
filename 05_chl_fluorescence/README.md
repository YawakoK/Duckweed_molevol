# 05_chl_fluorescence — Fv/Fm photoinhibition and recovery under high light

Produces: Fig. 4a, 4b; Supplementary Fig. 3; Supplementary Data 10–14
Author: Natsu Katayama

## Overview

| Step | Script | What it does | Output |
|---|---|---|---|
| 01 | `01_extract_fvfm.R` | raw PAM exports → HL (pfd = 1000) Fv/Fm at Day 0–4 for the six lineages; last 6 measurements per lineage × day | `output/SupplementaryData10_FvFm_extracted_data.csv` |
| 02 | `02_fvfm_stats.R` | Day 0/1/4: `lm(Fv_Fm ~ sample * time)`, `car::Anova(type = 2)`, emmeans Tukey time contrasts within lineage, Sidak compact letters among lineages within day | `output/SupplementaryData11_…` – `output/SupplementaryData14_…` |
| 03 | `03_fig4ab.R` | Fig. 4a (mean ± SE time course), Fig. 4b (Day 1 vs Day 4 boxplots, Tukey brackets) | `output/Fig4a.pdf`, `output/Fig4b.pdf` |
| 04 | `04_suppfig3.R` | Supplementary Fig. 3 (distributions by day with letters) and sessionInfo | `output/SupplementaryFig3.pdf`, `../environment/R_sessionInfo_05_chl_fluorescence.txt` |

| Output file | Comms Bio item |
|---|---|
| `SupplementaryData10_FvFm_extracted_data.csv` | Supplementary Data 10 (180 measurements, Day 0–4) |
| `SupplementaryData11_FvFm_summary_by_species_and_time.csv` | Supplementary Data 11 |
| `SupplementaryData12_FvFm_time_pairwise_comparisons_within_species.csv` | Supplementary Data 12 |
| `SupplementaryData13_FvFm_species_time_ANOVA.csv` | Supplementary Data 13 |
| `SupplementaryData14_FvFm_species_group_letters_within_time.csv` | Supplementary Data 14 |

The combined Fig. 4 consists of panels a and b from this folder and panels c and d from
`../06_pam_light_response`. Panels were assembled, and species names abbreviated and fonts
adjusted, in a vector editor; the R output here is the content of each panel before that step.

## Inputs (`data/`)

| File | Content | Source |
|---|---|---|
| `dat240408.txt`, `dat240520.txt`, `dat240527.txt`, `dat240610.txt`, `dat240617.txt`, `dat240624.txt`, `dat240722.txt`, `dat250420.txt`, `dat250609.txt` | PAM Fv/Fm exports (columns date, sample, F0, Fm, Fv_Fm, cal_Fv_Fm, Plate, hole, pfd, time, ld, medium); source of Sp7498, Lp9387, LaM10E, Wh9525, Wa | raw instrument tables (one per experiment start date) |
| `dat1.txt` – `dat5.txt` | earlier Fv/Fm time-course exports with time as elapsed days/hours (`1d`, `1d3h`, …); source of Lgp8L (*L. gibba*) | raw instrument tables |

Selection rules (step 01): pfd = 1000 (HL) only; Day 0–4; *L. gibba* is taken only from
`dat1`–`dat5` (elapsed time rounded to days, the 3 h point excluded); the last six
measurements in file order are kept per lineage × day (`slice_tail(n = 6)`). Other
accessions present in the raw tables (e.g. Lt6619, Si7178, Lm5512) are not used.

Lineage codes: Sp7498 *Spirodela polyrhiza*, Lp9387 *Landoltia punctata*, Lgp8L *Lemna gibba*,
LaM10E *Lemna aequinoctialis*, Wh9525 *Wolffiella hyalina*, Wa *Wolffia australiana*;
`antho` = anthocyanin-accumulating (Sp, Lp, Lg) vs `non_antho`.

## How to run

```bash
Rscript 01_extract_fvfm.R
Rscript 02_fvfm_stats.R
Rscript 03_fig4ab.R
Rscript 04_suppfig3.R
```

Runs in a few seconds. `set.seed(1)` fixes only the jitter of the plotted points in Fig. 4b and
Supplementary Fig. 3 (no random step in the analysis). The submitted figures had no fixed
seed, so point jitter differs slightly from the submitted PDFs.

## Values reported in the main text and where they are printed

| Main text | Value | Printed by |
|---|---|---|
| Methods "Statistical analysis" | Fv/Fm n = 180 measurements (Supplementary Data 10) | `01` |
| Results "Lineage-specific photoinhibition and recovery under high light" | lineage, time, lineage × time all P < 0.001 (F5,90 = 48.76, P = 3.6e-24; F2,90 = 569.0, P = 8.4e-52; F10,90 = 12.99, P = 9.5e-14) | `02` (Supplementary Data 13) |
| same | Day 1 → Day 4 increase significant in *S. polyrhiza* (P = 1.7e-8), *L. gibba* (P = 0.012), *W. australiana* (P = 0.0068); not in *L. punctata* (0.24), *L. aequinoctialis* (0.77), *W. hyalina* (0.93) | `02`, `03` (Supplementary Data 12, Fig. 4b) |
| same | Day 0 highest *W. australiana*, lowest *W. hyalina*; Day 1 highest *W. australiana*/*W. hyalina*, lowest *S. polyrhiza*/*L. gibba*; Day 4 *W. australiana* > *W. hyalina*, lowest *L. gibba* | `02` (Supplementary Data 11) |
| Supplementary Fig. 3 letters | compact letters within each day | `02` (Supplementary Data 14) |

## Notes

- The statistics use Day 0, 1 and 4 only (108 of the 180 measurements; residual df = 90).
  Day 2 and Day 3 are included in Supplementary Data 10 but not analysed.
- Multiple-comparison adjustment: the Day contrasts within lineage (Supplementary Data 12,
  Fig. 4b) are Tukey-adjusted. The compact letters (Supplementary Data 14, Supplementary Fig. 3)
  were computed with `multcomp::cld(..., adjust = "sidak")`, i.e. Sidak-adjusted pairwise
  comparisons, and the `lower.CL`/`upper.CL` columns of Supplementary Data 14 are
  Sidak-adjusted 95% CIs. With Tukey adjustment the letters are identical (checked in `02`).
- One Fv_Fm cell in `dat240520.txt` is recorded as `-`; it is coerced to NA (R warns) and is
  not among the selected measurements.

## Environment

Outputs were re-verified under R 4.6.0 (the manuscript reports R 4.4.1); dplyr 1.2.1,
readr 2.2.0, ggplot2 4.0.3, car 3.1-5, emmeans 2.0.3, multcomp 1.4-30.
Full list: `../environment/R_sessionInfo_05_chl_fluorescence.txt`.
