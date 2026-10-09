#!/bin/bash
# De novo transcriptome assembly with Trinity v2.15.2 (default settings: in-silico
# read normalisation to max. coverage 200, built-in salmon expression filter) for
# the three species without a reference genome.
#
# Read libraries per species (strand-specific, 150-bp paired-end):
#   Lemna aequinoctialis  9 libraries  (Azenta, NovaSeq 6000; BioProject PRJNA1508391)
#   Landoltia punctata    9 libraries  (Azenta, NovaSeq 6000; BioProject PRJNA1508391)
#   Wolffiella hyalina    9 libraries  (Azenta, this study)
#                       + 2 libraries sequenced in 2022 (this study, same BioProject)
# All libraries of a species were concatenated into left.fa / right.fa (FASTA) before
# assembly.  Trinity was run with the working directory on local disk (its
# Chrysalis/Butterfly stages create very many small files); salmon 1.10.0 was placed
# first on PATH because Trinity 2.15.2 is incompatible with the salmon 2.x index format.
#
# Usage: 01_trinity_assembly.sh <SPECIES> <dir with left.fa right.fa> <MAXMEM e.g. 200G> <CPU>
set -euo pipefail
SP=$1; IN=$2; MEM=$3; CPU=$4
Trinity --seqType fa --max_memory "$MEM" --CPU "$CPU" \
        --left "$IN/left.fa" --right "$IN/right.fa" \
        --output "${SP}_trinity"
# result: ${SP}_trinity.Trinity.fasta  (Trinity >= 2.15 writes the assembly next to the output dir)
