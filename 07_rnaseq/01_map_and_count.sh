#!/bin/bash
# Read mapping and counting for the S. polyrhiza / W. australiana high-light RNA-seq.
#
# 18 libraries per species: LL (WL_3h), HL_3h (SL_3h), HL_4d (SL_4D) x 6 replicates;
# replicates 1-3 sequenced by Azenta (Illumina NovaSeq 6000, files <sample>_R1/R2_001.fastq.gz)
# and replicates 4-6 by BGI (DNBSEQ-G400, files <sample>_1/2.fq.gz).  BioProject PRJNA1508391.
#
# Mapping references = nuclear genome + organellar genomes, with the corresponding
# combined GFF3 (nuclear + plastid + mitochondrial gene models, tRNA/rRNA removed);
# built in 01_sequence_prep/05_prepare_mapping_references.sh:
#   Sp: Sp9509_oxford_v3_withorga.fasta  + Sp_polyrhiza_9509-REF-OXFORD-3.0_CSHL2022v1.withorga_SP.gff3
#   Wa: Wo_australiana_8730-REF-CSHL-1.0_withorga.fasta + Wo_australiana_8730-REF-CSHL-1.0_CSHL2022v1.withorga_SP.gff3
#
# Tools as run: hisat2 2.2.1 (default parameters), featureCounts 2.0.8
#   featureCounts -p -t exon -g Parent -M --fraction -O
#   (-M --fraction: multi-mapping reads counted fractionally; -O: reads overlapping
#    several isoforms of one gene counted fractionally rather than discarded)
#
# Usage: 01_map_and_count.sh <Sp|Wa> <ref.fasta> <ref.gff3> <azenta_fastq_dir> <bgi_fastq_dir>
set -euo pipefail
SP=$1; REF=$2; GFF=$3; AZ=$4; BGI=$5
THREADS=${THREADS:-30}
SAM=sambam_$SP; mkdir -p "$SAM" output

[ -s "$REF.indexforhisat2.1.ht2" ] || hisat2-build "$REF" "$REF.indexforhisat2"

for cond in SL-3h SL-4D WL-3h; do
  for rep in 1 2 3; do                              # Azenta
    s=${SP}-${cond}-${rep}
    hisat2 -x "$REF.indexforhisat2" -1 "$AZ/${s}_R1_001.fastq.gz" -2 "$AZ/${s}_R2_001.fastq.gz" \
           -p "$THREADS" -S "$SAM/$s.sam"
  done
  c2=${cond//-/_}
  for rep in 1 2 3; do                              # BGI (samples 4-6 in the count matrix)
    s=${SP}_${c2}_${rep}
    hisat2 -x "$REF.indexforhisat2" -1 "$BGI/$s/${s}_1.fq.gz" -2 "$BGI/$s/${s}_2.fq.gz" \
           -p "$THREADS" -S "$SAM/$s.sam"
  done
done

# column order in the count matrix (= data/featureCounts_<SP>.txt):
#   Azenta SL-3h-1..3, SL-4D-1..3, WL-3h-1..3, then BGI SL_3h_1..3, SL_4D_1..3, WL_3h_1..3
featureCounts -p -T "$THREADS" -t exon -g Parent -M --fraction -O -a "$GFF" \
    -o "output/featureCounts_$SP.txt" \
    $SAM/${SP}-SL-3h-{1,2,3}.sam $SAM/${SP}-SL-4D-{1,2,3}.sam $SAM/${SP}-WL-3h-{1,2,3}.sam \
    $SAM/${SP}_SL_3h_{1,2,3}.sam $SAM/${SP}_SL_4D_{1,2,3}.sam $SAM/${SP}_WL_3h_{1,2,3}.sam
