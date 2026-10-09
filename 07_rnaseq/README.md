# 07_rnaseq — Lineage-specific gene expression in response to high light

Produces: Fig. 5; Supplementary Figs. 4–7; Supplementary Data 18, 19, 20
Author: Yawako W. Kawaguchi

## Overview

| Step | Script | What it does | Output |
|---|---|---|---|
| 01 | `01_map_and_count.sh` | hisat2 2.2.1 mapping to nuclear + organellar references; featureCounts 2.0.8 (`-p -t exon -g Parent -M --fraction -O`) | `data/featureCounts_Sp.txt`, `data/featureCounts_Wa.txt` |
| 02 | `02_build_orthogroup_maps.R` | orthogroup ↔ gene maps from `Orthogroups.tsv`; assigns the stable `OrthoID` | `output/Orth_{ath,At,Sp,Wa}.Rdata` |
| 03 | `03_tpm_batch_limma_anova.R` | TPM → sum per orthogroup → log2(TPM+1) → `limma::removeBatchEffect` (sequencing platform) → within-species limma contrasts → per-OG species × condition ANOVA (BH) | `output/logTPM_orthogroup_batch_corrected.*`, `output/limma_*.csv`, `output/anova_*.csv` |
| 04 | `04_transition_patterns.R` | U/F/D classification of LL→HL_3h and HL_3h→HL_4d (adj. *p* < 0.05, \|log2FC\| ≥ 0.5) → nine trajectories | `output/transition_pattern_pergene.csv`, `output/transition_pattern_9x9.csv` |
| 05 | `05_go_enrichment.R` | topGO 2.58.0 weight01 Fisher for 18 cohorts (species-biased, both-induced, trajectory) | `output/GO_*.csv` |
| 06 | `06_fig5.R` | Fig. 5 A–I | `output/Fig5.pdf` |
| 07 | `07_supfigs.R` | representative-orthogroup panels | `output/SupplementaryFig4_…pdf` – `output/SupplementaryFig7_…pdf` (6 = ROS-scavenging genes, 7 = *W. australiana*-biased) |
| 08 | `08_suptables.R` | supplementary data tables | `output/SupplementaryData20_GO_enrichment_all.csv`, `output/SupplementaryData19_per_orthogroup_response.csv`, `output/SupplementaryData18_orthogroups.csv` |

## Inputs (`data/`)

| File | Content | Source |
|---|---|---|
| `featureCounts_Sp.txt`, `featureCounts_Wa.txt` (+ `.summary`) | featureCounts matrices, 18 libraries each | step 01 (as run 2026-05-15) |
| `Orthogroups.tsv` | OrthoFinder v3.1.3 orthogroups, 8 species | `02_phylogenomics/01` |
| `Sp_gene_positions.txt`, `Wa_gene_positions.txt` | chr, start, end, gene from the combined GFF3s | `awk '$3=="gene"{print $1,$4,$5,$9}' <gff3> \| sed 's/ID=//'` |
| `ATH_GO_GOSLIM_gene2GO.txt` | *A. thaliana* gene → GO term (two columns) | TAIR `ATH_GO_GOSLIM.txt`, accessed 2026-02-23 (https://www.arabidopsis.org/download/list?dir=GO_and_PO_Annotations%2FGene_Ontology_Annotations), columns 1 and 5 |
| `Araport11_symbols.tsv` | Araport11 transcript → gene symbol | headers of `Araport11_pep_20250214_representative_gene_model` (TAIR) |

Sample naming inside the count matrices and all R objects:
`<Sp|Wa>_<WL_3h|SL_3h|SL_4D>_<1-6>` where **WL_3h = LL**, **SL_3h = HL_3h**, **SL_4D = HL_4d**;
replicates 1–3 were sequenced by Azenta (Illumina NovaSeq 6000) and 4–6 by BGI
(DNBSEQ-G400) — this is the `batch` factor removed in step 03.
Raw reads: NCBI BioProject PRJNA1508391.

## How to run

```bash
Rscript 02_build_orthogroup_maps.R        # ~10 s
Rscript 03_tpm_batch_limma_anova.R        # ~2 min (per-orthogroup ANOVA)
Rscript 04_transition_patterns.R          # ~5 s
Rscript 05_go_enrichment.R                # ~15 min (54 topGO runs)
Rscript 06_fig5.R                         # Fig. 5
Rscript 07_supfigs.R                      # Supplementary Figs. 4-7
Rscript 08_suptables.R                    # Supplementary Data 18, 19, 20
```

Step 01 needs the raw reads and the mapping references from `01_sequence_prep/05`;
its output is shipped in `data/`, so steps 02–08 run as is.

Values reported in the main text and where they are printed:

| Main text | Value | Printed by |
|---|---|---|
| L384 | 10,916 orthogroups expressed in both species | `03` (`orthogroups expressed in both species`), `04` (`N`) |
| L389–391 | Sp UF 225; Wa FU 655, FD 238, DF 245 | `04` (marginal pattern table); Fig. 5a |
| L391 | ~3.8 % identical trajectories among responsive OGs | `04` (`Identical EX FF/FF`) |
| Fig. 5 legend | FF class 95.8 % (Sp) / 88.2 % (Wa) | `04` |
| L397, L412, L420 | HL4d vs LL: 204 both-up, 1,792 Sp-biased, 435 Wa-biased | `06` (scatter counts), `08` (`Category_HL4d`) |
| L407–408 | HL3h vs LL: 167 Sp-biased, 29 Wa-biased (39 both-up) | `06` |
| L398–404, L409–431 | GO terms, FE and *p* | `output/SupplementaryData20_GO_enrichment_all.csv` (filter `comparison`, `cohort`) |
| Fig. 5d–i titles | species × condition interaction *p* (HL4d vs LL) | `06` |

## Environment

R 4.4.1; limma 3.62.2, topGO 2.58.0, GO.db 3.20.0, tidyverse 2.0.0, ggplot2 4.0.1,
patchwork 1.3.2, ggtext, ggh4x.  Full list: `../environment/R_sessionInfo_07_rnaseq.txt`.
