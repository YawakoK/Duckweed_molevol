# 03_gbif_distribution — Global distribution and latitude of anthocyanin-accumulating and non-accumulating duckweeds

Produces: Fig. 2a, 2b; Supplementary Fig. 1; Supplementary Data 1, 2
Author: Natsu Katayama

## Overview

| Step | Script | What it does | Output |
|---|---|---|---|
| 01 | `01_extract_columns.sh` | extracts 10 columns from the five GBIF SIMPLE_CSV downloads (awk) | `gbif_raw/<Genus>_extracted.txt` (not shipped) |
| 02 | `02_compact_occurrences.R` | combines the five extracts in a fixed order and keeps the 7 columns needed for filtering/thinning | `data/gbif_occurrences_compact.csv.xz` |
| 03 | `03_filter_and_thin.R` | species matching (36 accepted species), coordinate cleaning, duplicate removal, 0.5° spatial thinning with `set.seed(123)`; checks the result against Supplementary Data 2 | `output/filtering_counts.csv`, `output/thinned_records.csv` |
| 04 | `04_species_summary.R` | per-species latitude summary; median absolute latitude by anthocyanin status (all genera and *Lemna* only) | `output/SupplementaryData01_GBIF_species_and_latitude_summary.csv`, `output/abs_latitude_by_anthocyanin_status.csv` |
| 05 | `05_fig2.R` | Fig. 2a map, Fig. 2b absolute-latitude boxplot | `output/Fig2.pdf` |
| 06 | `06_suppfig1.R` | Supplementary Fig. 1 (*Lemna* only): map (page 1), absolute latitude (page 2); writes sessionInfo | `output/SupplementaryFig1.pdf` |

## Inputs (`data/`)

| File | Content | Source |
|---|---|---|
| `species_list_anthocyanin_status.tsv` | 36 accepted Lemnoideae species (`genus_cor`, `epithet`, `sp`) and anthocyanin status (`antho`: `Antho` = accumulating, `NonAntho` = non-accumulating; Les et al. 1997) | manual; two names carry a trailing space, removed in step 03 with `str_squish()` |
| `gbif_occurrences_compact.csv.xz` | all 753,829 raw GBIF records, in original order; columns `source_file`, `gbifID`, `species`, `scientificName`, `decimalLatitude`, `decimalLongitude`, `coordinateUncertaintyInMeters` | step 02 |
| `SupplementaryData02_GBIF_thinned_occurrence_records.csv` | **= Supplementary Data 2** (18,227 thinned records, identical file) | step 03 reproduces it exactly (see below) |

GBIF downloads (19 September 2024, SIMPLE_CSV; not shipped, ~380 MB in total).  Place them in
`gbif_raw/` named `<download key>.csv` to run steps 01–02:

| Taxon (GBIF taxon key) | Download key | DOI | Records |
|---|---|---|---|
| *Spirodela* (2867778) | 0022526-240906103802322 | https://doi.org/10.15468/dl.tqazwm | 125,396 |
| *Landoltia punctata* (2871369) | 0022527-240906103802322 | https://doi.org/10.15468/dl.cmkp53 | 3,960 |
| *Lemna* (2867567) | 0022523-240906103802322 | https://doi.org/10.15468/dl.sp7u55 | 601,484 |
| *Wolffiella* (2867437) | 0022531-240906103802322 | https://doi.org/10.15468/dl.wxhh5g | 2,669 |
| *Wolffia* (5330077) | 0022602-240906103802322 | https://doi.org/10.15468/dl.shuzby | 20,320 |

Note: *Landoltia* records also occur in the *Spirodela* download (as *Spirodela punctata*);
species are assigned from the first two words of `scientificName`, and the duplicate
species × coordinate records from the *Landoltia* download are removed in step 03.

## How to run

```bash
Rscript 03_filter_and_thin.R   # ~10 s
Rscript 04_species_summary.R
Rscript 05_fig2.R
Rscript 06_suppfig1.R
```

Steps 01–02 need the GBIF downloads (see above); their output is shipped in `data/`, so
steps 03–06 run as is.  Step 01 was re-run on the original downloads and reproduces the
extracts used for the paper byte-for-byte.

## Spatial thinning and reproducibility

Thinning keeps one randomly chosen record per species per 0.5° × 0.5° cell
(`round(coordinate * 2) / 2`, `dplyr::slice_sample(n = 1)` after `set.seed(123)`).
Step 03 reproduces Supplementary Data 2 exactly (same 18,227 `gbifID`s in the same row
order; checked under R 4.6.0, dplyr 1.2.1).  Because the random draw depends on the R / dplyr
sampling implementation and on the input row order, steps 04–06 read the shipped
Supplementary Data 2 (`data/`), not the step 03 output, so that the figures and statistics
are fixed to the submitted record set.  `output/thinned_records.csv` is kept only for the check.

Filtering counts printed by step 03 (`output/filtering_counts.csv`):

| Step | Records | Species |
|---|---|---|
| raw GBIF records | 753,829 | |
| species-level identification and coordinates present | 714,270 | |
| assigned to one of the 36 accepted species | 710,642 | 29 |
| valid coordinate range | 710,642 | 29 |
| coordinate uncertainty ≤ 10 km or missing | 703,837 | 29 |
| not at (0, 0) | 703,743 | 29 |
| species × coordinate duplicates removed | 443,314 | 29 |
| 0.5° spatial thinning | **18,227** | **29** |

Seven listed species have no GBIF records: *Lemna obscula*, *Wolffiella neotropica*,
*W. caudata*, *W. rotunda*, *W. brasiliensis*, *Wolffia neglecta*, *W. cylindraceae*.

## Values reported in the main text and where they are printed

| Main text | Value | Printed by |
|---|---|---|
| Results (Fig. 2) and Methods "GBIF" | 18,227 occurrence records from 29 species | `03` (last row of the counts table), `04` |
| Methods | 36 accepted duckweed species | `03` (`Accepted species in the list`) |
| Results (Fig. 2b) | median absolute latitude 48.3° (accumulating, 14,448 records / 9 species) vs 37.5° (non-accumulating, 3,779 / 20) | `04` (rows `All genera`) |
| Results, Discussion (Supplementary Fig. 1) | *Lemna*: 48.9° (accumulating, 10,591 / 6) vs 34.9° (non-accumulating, 2,212 / 5) | `04` (rows `Lemna`) |

## Verification against the submission

- Supplementary Data 2: reproduced exactly by step 03 (`gbifID`, species, coordinates, order).
- Fig. 2 and Supplementary Fig. 1: `output/Fig2.pdf` and `output/SupplementaryFig1.pdf` rasterise
  pixel-identically to the submitted PDFs (the PDF bytes differ only in metadata such as the creation date).
- Supplementary Data 1: `n_records` per species is identical to the submitted table, but
  `median_abs_latitude`, `mean_abs_latitude`, `min_latitude` and `max_latitude` differ slightly
  (e.g. *Landoltia punctata* median 33.598 vs 33.6).  The submitted Supplementary Data 1 was
  computed from a different random thinning draw than the one saved as Supplementary Data 2
  (the per-species record counts do not depend on the draw).  The version in `output/` is
  computed from Supplementary Data 2, consistent with Fig. 2 and the main-text medians.

## Environment

R 4.6.0 (the manuscript reports R 4.4.1; outputs were re-verified under R 4.6.0);
tidyverse 2.0.0 (dplyr 1.2.1, ggplot2 4.0.3, readr 2.2.0, stringr 1.6.0), maps, ggbeeswarm, patchwork.
Full list: `../environment/R_sessionInfo_03_gbif_distribution.txt`.
