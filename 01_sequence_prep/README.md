# 01_sequence_prep — Sequence data preparation

Produces: the protein and CDS sets used by `02_phylogenomics` (OrthoFinder input,
codon alignments) and the mapping references used by `07_rnaseq`; the de novo
assemblies and representative sequences deposited as supplementary sequence data.
Author: Yawako W. Kawaguchi

These are compute-heavy, one-off steps.  The scripts record the commands and
parameters as run; they are not expected to be re-executed by readers.  Their
products are available as follows.

| Product | Where |
|---|---|
| Trinity assemblies, representative protein/CDS (*L. aequinoctialis*, *La. punctata*, *W. hyalina*) | NCBI TSA under BioProject PRJNA1508391; also Zenodo dataset |
| *W. australiana* organellar gene models (MFannot mitochondrion + NCBI plastid) | Zenodo dataset |
| *S. polyrhiza*, *W. australiana*, *L. gibba* protein/CDS sets | derived from public Lemna Genome Hub assemblies with the scripts here |
| Raw RNA-seq reads (all species) | NCBI SRA, BioProject PRJNA1508391 |

## Scripts

| Script | What it does | Tools (version as run) |
|---|---|---|
| `01_trinity_assembly.sh` | de novo assembly per species | Trinity 2.15.2 (default; in-silico normalisation max cov 200), salmon 1.10.0 |
| `02_orf_prediction.sh` | CD-HIT-EST 98 % → TransDecoder.LongOrfs → Pfam-A (hmmscan) + Swiss-Prot (DIAMOND) → TransDecoder.Predict → representative isoform | CD-HIT 4.8.1, TransDecoder 5.5.0, HMMER 3.4, DIAMOND 2.2.5; Pfam-A (2023-12-26), Swiss-Prot (2023-12-27) |
| `03_select_representative_isoform.py` | one ORF per Trinity gene, highest coding score; protein and CDS kept in sync | Python 3.11 |
| `04_prepare_lemna_gibba_ref.sh` | AGAT standardisation and CDS/protein extraction from Le_gibba_7742a-REF-CSHL-1.0; longest transcript per gene (24,222) | AGAT 1.x |
| `05_prepare_mapping_references.sh` | organellar gene models (AGAT; MFannot for the *W. australiana* mitochondrion), combined nuclear + organellar genome FASTA / GFF3 for mapping, and the Sp / Wa protein and CDS sets | AGAT 1.0.0, gffread, seqkit, MFannot |

## Read libraries used for assembly

| Species | Libraries |
|---|---|
| *Lemna aequinoctialis* | 9 (Azenta, NovaSeq 6000, PE150; this study) |
| *Landoltia punctata* | 9 (Azenta, NovaSeq 6000, PE150; this study) |
| *Wolffiella hyalina* | 9 (Azenta, this study) + 2 libraries sequenced in 2022 (this study) |

Sequence counts of the resulting sets (= OrthoFinder input):

| Species | File | Sequences |
|---|---|---|
| *A. thaliana* | `Araport11_pep_20250214.fasta` | 48,354 peptides (all isoforms of 27,650 loci; Araport11 release 2025-02-14) |
| *C. esculenta* | `Colocasia_esculenta.faa` | 28,695 (distributed protein set, CNGBdb CNP0001082) |
| *S. polyrhiza* | `Sp_withorga.fasta` | 20,464 (primary nuclear proteins + plastid + mitochondrion) |
| *La. punctata* | `La_longest.isoform.faa` | 26,873 |
| *L. aequinoctialis* | `Laeq_longest.isoform.faa` | 18,630 |
| *L. gibba* | `Lgib.faa` | 24,222 |
| *W. hyalina* | `Wo_longest.isoform.faa` | 33,332 |
| *W. australiana* | `Wa_withorga.fasta` | 21,022 (primary nuclear proteins + plastid + MFannot mitochondrion) |

Label caveat used throughout the repository: `La` = *Landoltia punctata*,
`Laeq` = *Lemna aequinoctialis*, `Wo`/`Wh` = *Wolffiella hyalina*, `Wa` = *Wolffia australiana*.
