#!/bin/bash
# Orthogroup inference across eight species with OrthoFinder v3.1.3 (default
# parameters).  Produces Orthogroups.tsv / Orthogroups.GeneCount.tsv used by every
# downstream step of this directory and by 07_rnaseq.
#
# Input directory: one protein FASTA per species (one representative sequence
# per gene).  File names become the species column names in Orthogroups.tsv:
#
#   Araport11_pep_20250214.fasta   Arabidopsis thaliana  (Araport11, TAIR)
#   Colocasia_esculenta.faa        Colocasia esculenta   (CNGBdb CNP0001082)
#   Sp_withorga.fasta              Spirodela polyrhiza   (9509-REF-OXFORD-3.0 + organelles; 01_sequence_prep/05)
#   La_longest.isoform.faa         Landoltia punctata    (de novo transcriptome; 01_sequence_prep/02)
#   Laeq_longest.isoform.faa       Lemna aequinoctialis  (de novo transcriptome; 01_sequence_prep/02)
#   Lgib.faa                       Lemna gibba           (7742a-REF-CSHL-1.0; 01_sequence_prep/04)
#   Wo_longest.isoform.faa         Wolffiella hyalina    (de novo transcriptome; 01_sequence_prep/02)
#   Wa_withorga.fasta              Wolffia australiana   (8730-REF-CSHL-1.0 + organelles; 01_sequence_prep/05)
#
# As run (2026-08-24):
#   orthofinder -f Eight_sp -t 28 -a 8 -o orthofinder3_output
# Search: diamond; MSA: famsa; tree: fasttree (OrthoFinder 3 defaults).
set -euo pipefail
IN=${1:-Eight_sp}
OUT=${2:-orthofinder3_output}
orthofinder -f "$IN" -t "${THREADS:-28}" -a "${THREADS_ALG:-8}" -o "$OUT"
