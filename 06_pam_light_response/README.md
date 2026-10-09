# 06_pam_light_response — PAM light-response curves and PSII energy partitioning

Produces: Fig. 4c, Fig. 4d; Supplementary Data 15, 16, 17
Author: Natsu Katayama

Fig. 4 in the paper combines panels a and b (Fv/Fm; folder `05_chl_fluorescence`)
with panels c and d from this folder. The final panel arrangement, species
abbreviations (*Sp*, *Wa*) and P-value typography were set in the layout software
(Affinity); the data, colours, shapes and line types are those produced here.

## Overview

| Step | Script | What it does | Output |
|---|---|---|---|
| 01 | `01_process_light_curves.R` | reads the 12 raw light-curve exports (*S. polyrhiza* Sp7498 and *W. australiana* Wa8730, LL and HL, 3 plant samples each), attaches the measured PPFD of each of the 17 actinic light steps, computes ETR = Y(II) × PAR × 0.84 × 0.5 | `output/PAM_light_curves_processed.csv` |
| 02 | `02_quantum_yield_stats.R` | extracts Y(II), Y(NPQ), Y(NO) at AL50 and AL1000; summary; LMM `Y ~ species * growth_light * AL + (1 \| sample)` per yield with Sidak-adjusted Sp − Wa contrasts of estimated marginal means (emmeans); per-AL `lm(Y ~ species * growth_light)` interaction ANOVA | `output/SupplementaryData15_PAM_quantum_yield_summary.csv`, `output/SupplementaryData16_PAM_species_pairwise_contrasts_by_condition.csv`, `output/SupplementaryData17_PAM_species_growth_light_interaction_ANOVA.csv` |
| 03 | `03_fig4c_light_response_curves.R` | Fig. 4c: mean Y(II), ETR, NPQ vs actinic PPFD | `output/Fig4c.pdf` |
| 04 | `04_fig4d_quantum_yield_partitioning.R` | Fig. 4d: mean ± SE of the three yields at AL50/AL1000, LL vs HL, with Sp vs Wa brackets (P from Supplementary Data 16); writes sessionInfo | `output/Fig4d.pdf` |

## Inputs (`data/`)

| File | Content | Source |
|---|---|---|
| `raw_PAM_light_curves/{Sp,Wa}_50_rep{1,2,3}_*.CSV` | raw light-curve exports, LL-grown plants (growth PPFD 50), measured 2024-07-02 to 07-04 | PAM fluorometer export (`;`-separated) |
| `raw_PAM_light_curves/{Sp,Wa}_1000_DAY4_rep{1,2}_*.CSV`, `{Sp,Wa}_1000_LC_rep3_*.CSV` | raw light-curve exports, HL-grown plants (growth PPFD 1000, day 4), measured 2024-07-12 and 2024-07-25 | as above |
| `raw_PAM_light_curves/PARlist_240702.csv` | measured PPFD (µmol m⁻² s⁻¹) of the 17 actinic steps of the light-curve program | measured for the light-curve program, 2024-07-02 |

Light-curve measurements of four other lineages made in the same campaign, and an
additional LL replicate 4 of Sp and Wa (2024-07-05; excluded because the growth-light
intensity on that day was suspected to be higher than the LL setting), are not used in
the paper and are not shipped.

Definitions used in the scripts:

- growth light `PFD_50` = **LL**, `PFD_1000` = **HL**;
- **AL50** = light step 2, measured PPFD 61.885; **AL1000** = light step 12, measured PPFD 1012.5
  (internal labels `PPFD_50`, `PPFD_1000`);
- ETR uses the actinic PAR value recorded by the instrument for each step (column `PAR`
  of the raw export); the x axis of Fig. 4c uses the measured PPFD from `PARlist_240702.csv`.
- **Sample size**: n = 3 plant samples (independent light curves) per species × growth
  light × AL combination (12 light curves in total). This is the value missing in the
  manuscript's Statistical analyses ("n = ?"); it is printed by step 02 and is column `n`
  of Supplementary Data 15.

## How to run

```bash
Rscript 01_process_light_curves.R
Rscript 02_quantum_yield_stats.R
Rscript 03_fig4c_light_response_curves.R
Rscript 04_fig4d_quantum_yield_partitioning.R
```

No random numbers are used. All steps run in a few seconds.

## Values reported in the main text and where they are printed

| Main text | Value | Printed by |
|---|---|---|
| Results "Photosynthetic light-response curves and PSII energy partitioning" | Y(II) Wa > Sp: AL50 LL P = 0.005; AL50 HL P < 0.001; AL1000 LL and HL P < 0.001 | `02` (Supplementary Data 16 table) |
| same paragraph | Y(NPQ) Sp > Wa at AL1000, HL-grown, P = 0.002 | `02` (Supplementary Data 16 table) |
| same paragraph | no species × growth-light interaction for any yield at either AL (all P > 0.05; smallest P = 0.081, Y(II) at AL1000) | `02` (Supplementary Data 17 table, `min interaction P`) |
| Methods "Statistical analyses" | n = 3 per species × growth light × AL | `02` (`n per species x growth light x AL`) |

## Verification

`output/SupplementaryData15/16/17_*.csv` were regenerated from the raw exports and
compared with the submitted Supplementary Data 15–17: identical columns, rows and labels;
numeric differences ≤ 1e-9 (floating-point only). The Sp/Wa rows of the processed
light-curve table are identical to the table used for the submitted figure.
Outputs were re-verified under R 4.6.0 (the manuscript states R 4.4.1).

## Environment

R 4.6.0; lme4 2.0.1, lmerTest 3.2.1, emmeans 2.0.3, pbkrtest 0.5.5 (Kenward–Roger d.f.,
emmeans default for `lmer` models), dplyr 1.2.1, tidyr 1.3.2, ggplot2 4.0.3.
Full list: `../environment/R_sessionInfo_06_pam_light_response.txt`.
