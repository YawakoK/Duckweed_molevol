# Shifts in high-light acclimation strategies associated with molecular evolutionary rates in duckweeds

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23253799.svg)](https://doi.org/10.5281/zenodo.23253799)

Analysis code and derived data for

> Kawaguchi Y.W., Kono M., Isoda M., Kawade K., Tabeta H., Sasaki R., Oikawa A., Hirai M.Y.,
> Josephine E., Valerie L., Katayama N.  *Shifts in high-light acclimation strategies associated
> with molecular evolutionary rates in duckweeds.* 

Each numbered directory corresponds to one results section / figure and is self-contained:
`README.md` (figure → script → output map), numbered scripts in execution order, `data/`
(inputs, ≤ 10 MB each) and `output/` (what the scripts produce, as submitted).

| Directory | Figures / tables | Analysis | Owner |
|---|---|---|---|
| `01_sequence_prep/` | — (sequence data) | Trinity assembly, ORF prediction, reference preparation | Y.W.K. |
| `02_phylogenomics/` | Fig. 1b–e | OrthoFinder, codon alignments, gene / species trees, root-to-tip rates, dN/dS (HyPhy) | Y.W.K. |
| `03_gbif_distribution/` | Fig. 2; Supplementary Fig. 1; Supplementary Data 1–2 | GBIF occurrence filtering, spatial thinning, latitude | N.K. |
| `04_growth_rate/` | Fig. 3; Supplementary Fig. 2; Supplementary Data 3–9 | frond area (OpenCV), RGR, growth-curve models, Gamma GLMM | N.K. |
| `05_chl_fluorescence/` | Fig. 4a–b; Supplementary Fig. 3; Supplementary Data 10–14 | Fv/Fm linear model, Tukey / Sidak contrasts | N.K. |
| `06_pam_light_response/` | Fig. 4c–d; Supplementary Data 15–17 | light-response curves, Y(II)/Y(NPQ)/Y(NO) mixed models | N.K. |
| `07_rnaseq/` | Fig. 5; Supplementary Figs. 4–7; Supplementary Data 18–20 | hisat2/featureCounts, orthogroup TPM, limma, ANOVA, trajectories, topGO | Y.W.K. |
| `08_metabolome/` | Fig. 6; Supplementary Fig. 8; Supplementary Data 21–22 | metabolite log2FC, lineage × light models | N.K. |

Supplementary Data 23 (sequence-data and RNA-seq sample information) is a metadata table and has no generating script.

## Quick start

Figures that depend only on the tables shipped in `data/` regenerate in minutes:

```bash
cd 02_phylogenomics && Rscript 08_fig1c_relative_rate.R && python3 09_fig1de_dnds.py   # Fig. 1c-e
cd ../07_rnaseq && for s in 02 03 04 05 06 07 08; do Rscript ${s}_*.R; done             # Fig. 5, Supplementary Figs. 4-7, Supplementary Data 18-20
cd ../03_gbif_distribution && for s in 03 04 05 06; do Rscript ${s}_*.R; done           # Fig. 2, Supplementary Fig. 1, Supplementary Data 1-2
cd ../04_growth_rate && for s in 03 04 05 06 07; do Rscript ${s}_*.R; done              # Fig. 3, Supplementary Fig. 2, Supplementary Data 4-9
cd ../05_chl_fluorescence && for s in 01 02 03 04; do Rscript ${s}_*.R; done            # Fig. 4a-b, Supplementary Fig. 3, Supplementary Data 10-14
cd ../06_pam_light_response && for s in 01 02 03 04; do Rscript ${s}_*.R; done          # Fig. 4c-d, Supplementary Data 15-17
cd ../08_metabolome && for s in 01 02 03 04 05; do Rscript ${s}_*.R; done               # Fig. 6, Supplementary Fig. 8, Supplementary Data 21-22
```

Every script is run from its own directory and uses relative paths only.  Where a
value appears in the main text, the script that prints it is listed in the
directory README.

## Data

| Data | Repository |
|---|---|
| Raw RNA-seq reads, transcriptome assemblies | NCBI BioProject [PRJNA1508391](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1508391) |
| Per-orthogroup codon alignments, gene trees, aBSREL JSONs (~1.4 GB); *W. australiana* organellar gene models | Zenodo dataset [10.5281/zenodo.23254999](https://doi.org/10.5281/zenodo.23254999) |
| GBIF occurrence downloads (19 September 2024) | https://doi.org/10.15468/dl.sp7u55 (*Lemna*), https://doi.org/10.15468/dl.tqazwm (*Spirodela*), https://doi.org/10.15468/dl.cmkp53 (*Landoltia*), https://doi.org/10.15468/dl.wxhh5g (*Wolffiella*), https://doi.org/10.15468/dl.shuzby (*Wolffia*) |
| Public genomes | Lemna Genome Hub (Sp 9509-REF-OXFORD-3.0, Wa 8730-REF-CSHL-1.0, Lg 7742a-REF-CSHL-1.0); CNGBdb CNP0001082 (*C. esculenta*); TAIR Araport11 |

## Environment

`environment/` holds the R `sessionInfo()` of each analysis (`R_sessionInfo_<dir>.txt`), the Python package list
(`python_requirements.txt`) and the versions of every command-line tool (`tools.md`).

## License

Code: MIT (`LICENSE`).  Data tables and figures in `data/` and `output/`: CC BY 4.0.

## Citation

See `CITATION.cff`.  Archived versions of this repository are available at Zenodo:
concept DOI [10.5281/zenodo.23253799](https://doi.org/10.5281/zenodo.23253799) (always resolves to the latest version);
v1.0.0 (initial submission) [10.5281/zenodo.23253800](https://doi.org/10.5281/zenodo.23253800).
