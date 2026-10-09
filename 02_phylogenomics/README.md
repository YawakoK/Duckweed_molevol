# 02_phylogenomics — Phylogeny and lineage-specific molecular evolutionary rates

Produces: Fig. 1b, 1C, 1D, 1E; Supplementary Data 18 (orthogroup membership, see 07_rnaseq/08)
Author: Yawako W. Kawaguchi

## Overview

| Step | Script | What it does | Output |
|---|---|---|---|
| 01 | `01_orthofinder.sh` | OrthoFinder v3.1.3 on eight proteomes | `Orthogroups.tsv`, `Orthogroups.GeneCount.tsv` (→ `data/`) |
| 02 | `02_build_orthogroup_cds.py` | select orthogroups (no duplicates, ≥ 6 of 7 taxa) and write per-OG CDS; flag set A (all 7 single copy, n = 2,819) | `TreeEachOG/<OG>/all.cds.fa`, `gene_sets.json` |
| 03 | `03_align_and_tree_one_og.sh` | per OG: translate → MAFFT → back-translate → IQ-TREE 2 gene tree (TIM3+F+I, 1000 UFBoot, *Colocasia* outgroup) | `TreeEachOG/<OG>/all.cds.aln.fa`, `<OG>.treefile` |
| 04 | `04_species_tree.sh` | AMAS concatenation of the 2,819 single-copy codon alignments → IQ-TREE 3 partitioned ML (1000 UFBoot) | `data/species_tree.treefile` (**Fig. 1b**, drawn with PhyloWeaver) |
| 05 | `05_absrel_one_og.sh` | HyPhy 2.5.61 aBSREL per OG (MG94×REV baseline + adaptive branch-site fit), stop codons masked | `absrel/absrel_json/<OG>.json` |
| 06 | `06_summarise_absrel.py` | per-branch table with gene-level QC; root-to-tip relative rates from the gene trees | `absrel_per_branch.csv`, `fig1_relative_rates.csv` |
| 07 | `07_roottotip_dnds.py` | root-to-tip dS and dN from the duckweed MRCA | `roottotip_dnds.csv` |
| 08 | `08_fig1c_relative_rate.R` | Friedman + pairwise Wilcoxon (Holm), compact letters, violins | `output/Fig1C_relative_rate.pdf`, `output/Fig1C_medians_letters.csv` (**Fig. 1c**) |
| 09 | `09_fig1de_dnds.py` | QC ∩ single-copy (n = 2,474); root-to-tip dS violins; median dN vs dS with proportional fit | `output/Fig1DE_dN_dS.pdf` (**Fig. 1d, 1E**) |

`run_pipeline.sh` chains steps 02 → 03 → 05 → 06 → 07.

## Inputs

* Proteomes for step 01 and per-species CDS for step 02 come from `01_sequence_prep/`
  (three de novo transcriptomes, *L. gibba* REF-CSHL-1.0 via AGAT, *S. polyrhiza* and
  *W. australiana* nuclear + organellar models) plus the public *A. thaliana* (Araport11
  2025-02-14) and *C. esculenta* (CNGBdb CNP0001082) proteomes.
* `data/Orthogroups.tsv`, `data/Orthogroups.GeneCount.tsv` — OrthoFinder output (step 01, as run 2026-08-24).
* `data/gene_sets.json` — orthogroup sets chosen in step 02 (`selected` n = 3,615; `setA` n = 2,819).
* `data/species_tree.treefile`, `data/species_tree.iqtree`, `data/species_tree_partitions.txt` — step 04 output.

## Checkpoint tables shipped in `data/`

Steps 03 and 05 take days of CPU and their outputs (3,615 codon alignments, gene trees and
aBSREL JSONs, ~1.4 GB) are archived separately (Zenodo dataset, see the top-level README).
Their summaries are shipped here so that the figures can be regenerated directly:

| File | Made by | Rows |
|---|---|---|
| `data/absrel_per_branch.csv` | 06 | 38,168 branches of 3,615 OGs; `passed_QC` = 1 for 3,116 OGs |
| `data/fig1_relative_rates.csv` | 06 | 2,819 set-A OGs × 6 species |
| `data/roottotip_dnds.csv` | 07 | 18,171 species × OG rows (QC-passing OGs) |

The figure scripts read these from `data/` by default; after re-running steps 06–07,
`IN=output Rscript 08_fig1c_relative_rate.R` / `IN=output python3 09_fig1de_dnds.py`
use the regenerated copies instead.

## How to run (figures only)

```bash
Rscript 08_fig1c_relative_rate.R      # Fig. 1c  (~5 s)
python3 09_fig1de_dnds.py             # Fig. 1d, 1E  (~5 s)
```

Values reported in the main text and where they are printed:

| Main text | Value | Printed by |
|---|---|---|
| L159, L195 | 2,819 single-copy orthogroups | `06` (`gene trees used`), `data/gene_sets.json` `setA` |
| L170 | all species pairs differ, *p* < 0.05 (Friedman + Wilcoxon/Holm) | `08` |
| L172–173 | median relative branch length 0.250 (*S. polyrhiza*) vs 1.260 (*W. hyalina*) | `08`, `output/Fig1C_medians_letters.csv` |
| L176 | 2,474 orthogroups passing QC | `09` (`orthogroups used`) |
| L181 | proportional fit ω = 0.50 | `09` |
| L182–183 | per-lineage ω 0.43 (*L. gibba*) – 0.56 (*S. polyrhiza*) | `09` |

## How to run (full pipeline)

Requires OrthoFinder 3.1.3, MAFFT 7.520, IQ-TREE 2.2.6 (gene trees) and 3.0.1
(species tree), HyPhy 2.5.61, AMAS, Python 3.11.

```bash
bash 01_orthofinder.sh Eight_sp orthofinder3_output
OGDIR=orthofinder3_output/Results_*/Orthogroups CDSDIR=<per-species CDS> NPROC=26 bash run_pipeline.sh
TREES=TreeEachOG bash 04_species_tree.sh
```

Tip labels used throughout: `Sp` *Spirodela polyrhiza*, `La` *Landoltia punctata*,
`Lgib` *Lemna gibba*, `Laeq` *Lemna aequinoctialis*, `Wa` *Wolffia australiana*,
`Wo` *Wolffiella hyalina*, `Colo` *Colocasia esculenta* (outgroup).  Note that `La`
is *Landoltia*, not *Lemna aequinoctialis* (`Laeq`), and `Wo` is *Wolffiella*, not
*Wolffia* (`Wa`).  The figure scripts relabel to the manuscript abbreviations
(Sp, Lp, Lg, La, Wa, Wh).

## Environment

Python 3.11.5; matplotlib 3.7.2, numpy 1.24.3 (step 09).
R 4.4.1; tidyverse 2.0.0 (step 08); full list in `../environment/R_sessionInfo_02_phylogenomics.txt`.
Command-line tool versions: `../environment/tools.md`.
