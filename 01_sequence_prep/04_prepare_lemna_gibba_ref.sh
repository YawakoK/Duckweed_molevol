#!/bin/bash
# Lemna gibba protein and CDS sets from the chromosome-level assembly
# Le_gibba_7742a-REF-CSHL-1.0 (Lemna Genome Hub; Ernst et al. 2025).
#
#   1. AGAT agat_convert_sp_gxf2gxf.pl   standardise the MAKER gff3
#   2. AGAT agat_sp_extract_sequences.pl  CDS (-t cds) and protein (-t cds -p)
#   3. keep the LONGEST transcript per gene (transcript IDs are <gene>_T<nnn>)
#
# The same AGAT extraction (steps 1-2) was applied to the S. polyrhiza and
# W. australiana annotations in 05_prepare_mapping_references.sh.
#
# Usage: 04_prepare_lemna_gibba_ref.sh <genome.fasta(.gz)> <genes.gff3>
# Output: Lgib.faa (OrthoFinder input), Lgib.cds.fa (codon analyses), 24,222 genes each
# Version as run: AGAT 1.x (conda), 2026-07-22.
set -euo pipefail
GENOME=$1; GFF=$2

echo "=== 1. genome ==="
case "$GENOME" in *.gz) zcat "$GENOME" > genome.fa ;; *) cp "$GENOME" genome.fa ;; esac

echo "=== 2. AGAT standardise ==="
agat_convert_sp_gxf2gxf.pl -g "$GFF" -o std.gff3 > agat_convert.log 2>&1

echo "=== 3. AGAT extract CDS and protein ==="
agat_sp_extract_sequences.pl -g std.gff3 -f genome.fa -t cds    -o all.cds.fa > agat_cds.log 2>&1
agat_sp_extract_sequences.pl -g std.gff3 -f genome.fa -t cds -p -o all.pep.fa > agat_pep.log 2>&1

echo "=== 4. longest isoform per gene ==="
python3 - <<'PY'
import re
GENE = re.compile(r"^(.*)_T\d+$")
def read(path):
    n = None; buf = []
    for line in open(path):
        if line.startswith(">"):
            if n is not None: yield n, "".join(buf)
            n = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if n is not None: yield n, "".join(buf)
cds = dict(read("all.cds.fa"))
best = {}
for sid, s in cds.items():
    m = GENE.match(sid); g = m.group(1) if m else sid
    if g not in best or len(s) > len(cds[best[g]]): best[g] = sid
keep = {best[g] for g in best}
print(f"  transcripts {len(cds)} -> genes {len(best)}")
for src, dst in (("all.pep.fa", "Lgib.faa"), ("all.cds.fa", "Lgib.cds.fa")):
    n = 0
    with open(dst, "w") as out:
        for sid, s in read(src):
            if sid in keep:
                out.write(f">{sid}\n{s}\n"); n += 1
    print(f"  {dst}: {n}")
PY
