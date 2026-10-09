# 08_metabolome — Lineage-specific metabolite responses to high light

Produces: Fig. 6a–e; Supplementary Fig. 8; Supplementary Data 21, 22
Author: Natsu Katayama

## Overview

| Step | Script | What it does | Output |
|---|---|---|---|
| 01 | `01_read_abundance.R` | reads the three per-sample abundance matrices (CE-TOF MS, GC-MS/MS, LC-MS/MS), removes unidentified features, sets missing CE/GC values to 0, reshapes to long format | `output/metabolite_abundance_long.csv` |
| 02 | `02_pca.R` | PCA of log2(value + 1), centred and scaled, per platform | `output/SupplementaryFig8.pdf`, `output/SupplementaryFig8_PCA_scores.csv` |
| 03 | `03_log2fc_interaction.R` | metabolite-specific pseudocount (half the smallest positive value) → log2(value + pseudocount) → log2FC = mean log2 HL − mean log2 LL per species → `lm(log2y ~ Lineage * Light)`, interaction F test → BH FDR within platform → `Sp > Wa` / `Wa > Sp` classes | `output/SupplementaryData21_all_metabolite_log2FC_statistics.csv`, `output/metabolite_log2_abundance_long.csv` |
| 04 | `04_categories_within_species.R` | name-based functional categories; one-sided Welch t-test HL > LL within each species (BH within platform and species); table of the 115 metabolites in the Fig. 6 categories | `output/SupplementaryData22_category_metabolite_log2FC_statistics.csv`, `output/metabolite_category_annotation.csv` |
| 05 | `05_fig6.R` | Fig. 6a–e: category log2FC scatter (Sp vs Wa) + log2 abundance of three representative metabolites per category | `output/Fig6.pdf`, `output/Fig6_representative_metabolites_source_data.csv` |

`output/Fig6.pdf` holds the R-generated panels; in the submitted Fig. 6 the same
panels were re-arranged (scatter on top, abundance plots below, one column per
category) in Affinity Designer without changing the data.

## Inputs (`data/`)

| File | Content | Source |
|---|---|---|
| `CEMS_abundance.csv` | CE-TOF MS relative abundance, 259 identified features × 20 samples; columns: metabolite (row name, `T-xxx_<name>`), `Ion` (Cation/Anion), samples, `MS` | CE-TOF MS quantification table, after removing "Unknown" features and setting NA to 0 |
| `GCMS_abundance.csv` | GC-MS/MS relative abundance, 167 identified features × 20 samples; columns: metabolite (row name, TMS-derivative name), samples, `MS` | GC-MS/MS quantification table, after removing "Unknown" features and setting NA to 0 |
| `LCMS_abundance.csv` | LC-MS/MS preprocessed peak areas, 1,304 identified features; BGI annotation columns (ion mode, RT, m/z, formula, `Name`, identification level, KEGG/HMDB IDs, classes, pathways), 20 samples and 9 pooled QC injections (QC not used) | BGI Genomics `2.1.DataPreprocessed` table, rows with `Name == "unidentified"` removed |

Sample columns are `<Sp|Wa>_<HL|LL>_rep<1-5>` (five biological replicates per
species × light; Day 4 samples). In the BGI table the original sample names
were `SpSL1–5` (Sp HL), `SpWL1–5` (Sp LL), `WaSL1–5` (Wa HL), `WaWL1–5` (Wa LL).

These three files are the processed matrices used for all analyses; they were
written from the platform quantification tables by removing unidentified features
(and, for CE/GC, replacing NA by 0). Step 01 re-applies the same rule. Raw instrument data are not part of this repository;
they will be deposited in a public metabolomics repository (repository TBD, e.g.
MetaboLights). LC-MS/MS metabolite extraction, peak extraction, identification,
preprocessing and QC were done by BGI Genomics; GC-MS/MS and CE-TOF MS processing
followed refs. 81–84 of the manuscript.

## How to run

```bash
cd 08_metabolome
Rscript 01_read_abundance.R              # ~2 s
Rscript 02_pca.R                         # Supplementary Fig. 8
Rscript 03_log2fc_interaction.R          # Supplementary Data 21 (~5 s, 1,730 linear models)
Rscript 04_categories_within_species.R   # Supplementary Data 22
Rscript 05_fig6.R                        # Fig. 6
```

No random numbers are used in the statistics. `05_fig6.R` fixes the horizontal
jitter of replicate points (`seed = 1`) and the label placement (`seed = 7`).

Values reported in the main text and where they are printed:

| Main text (Results "Lineage-specific metabolite accumulation in response to high light") | Value | Printed by |
|---|---|---|
| metabolites analysed (Supplementary Data 21) | 259 CE-TOF MS + 167 GC-MS/MS + 1,304 LC-MS/MS = 1,730 | `01`, `03` |
| flavonoid-related metabolites in Fig. 6a | 61 | `04` (`Metabolites per Fig. 6 panel`) |
| flavonoids with larger HL response in Sp / Wa | 18 / 12 | `04` (panel × `plot_class` table) |
| flavonoids significantly increased under HL in both species (all flavonol/flavone-related) | 6 | `04` |
| species-biased flavonoids also increased within species | 16 of 18 (Sp), 11 of 12 (Wa) | `04` |
| amino acids/peptides Wa > Sp / Sp > Wa | 14 / 1 | `04` |
| central carbon, sugars/carbohydrates, energy/redox: Wa > Sp (Sp > Wa) | 6, 10, 4 (0, 0, 0) | `04` |
| Wa-biased primary metabolites also increased within Wa | 30 of 34 | `04` |
| single Sp-biased primary metabolite also increased within Sp | see note below | `04` |
| interaction FDR of the 15 representative metabolites (Fig. 6a–e labels) | | `05` |
| PCA % variance (Supplementary Fig. 8 axes) | | `02` |

Interaction classes over all metabolites (FDR < 0.05; `03`): CE-TOF MS 60 Sp > Wa /
56 Wa > Sp; GC-MS/MS 10 / 83; LC-MS/MS 289 / 196.

## Notes

- `output/SupplementaryData21_*.csv` and `output/SupplementaryData22_*.csv` reproduce
  Supplementary Data 21 and 22 (same rows, columns and order; numeric differences
  < 1e-9 relative, from CSV round-tripping in the original workflow).
- FDR (interaction) and the within-species FDRs are computed across all metabolites
  of a platform, not only within the Fig. 6 categories.
- Categories are assigned by regular expressions on the metabolite names
  (`04_categories_within_species.R`); Fig. 6d combines the GC-MS/MS categories
  "Carbohydrates / Sugar metabolism" and "Polyols". The category of every metabolite
  is listed in `output/metabolite_category_annotation.csv`.
- The single Sp > Wa metabolite in Fig. 6b–e is asparagine (log2FC Sp −0.31,
  Wa −2.35); it decreased less in Sp than in Wa and is not significantly increased
  within Sp (one-sided FDR ≈ 1). `04` prints `Panels b-e, Sp > Wa also increased
  within Sp: 0 of 1`.
- Outputs were re-verified under R 4.6.0 (the manuscript states R 4.4.1).

## Environment

R 4.6.0; dplyr, tidyr, readr, stringr (tidyverse 2.0.0), ggplot2, patchwork 1.3.2,
ggrepel 0.9.8. Full list: `../environment/R_sessionInfo_08_metabolome.txt`.
