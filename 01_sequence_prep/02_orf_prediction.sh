#!/bin/bash
# Post-assembly ORF pipeline for one de novo transcriptome (Laeq | La | Wo).
#
#   1. CD-HIT-EST v4.8.1        -c 0.98 -p 1 -d 0 -b 3       (cluster at 98 % identity)
#   2. TransDecoder.LongOrfs v5.5.0
#   3. hmmscan (HMMER 3.4) vs Pfam-A         -E 1e-10 --domE 1e-10
#   4. DIAMOND v2.2.5 blastp vs Swiss-Prot   --max-target-seqs 1 --evalue 1e-10
#   5. TransDecoder.Predict   (retain ORFs supported by Pfam or Swiss-Prot homology)
#   6. one representative isoform per Trinity gene, highest coding score
#      (03_select_representative_isoform.py; protein and CDS written together)
#
# Usage: 02_orf_prediction.sh <SPECIES> <assembly.fasta> <Pfam-A.hmm> <swissprot.dmnd> [threads]
# Output: <SPECIES>_longest.isoform.faa (OrthoFinder input), <SPECIES>.cds.fa (codon analyses)
set -euo pipefail
SP=$1; ASM=$2; PFAM=$3; SPROT=$4; T=${5:-8}
SELECT=$(dirname "$0")/03_select_representative_isoform.py

echo "=== [$SP] 1. CD-HIT-EST ==="
cd-hit-est -i "$ASM" -o ${SP}_cdhit.fasta -c 0.98 -p 1 -d 0 -b 3 -T "$T" -M 16000

echo "=== [$SP] 2. TransDecoder.LongOrfs ==="
TransDecoder.LongOrfs -t ${SP}_cdhit.fasta

echo "=== [$SP] 3. HMMER vs Pfam-A ==="
hmmscan --cpu "$T" -E 1e-10 --domE 1e-10 --domtblout pfam.domtblout "$PFAM" \
        ${SP}_cdhit.fasta.transdecoder_dir/longest_orfs.pep > /dev/null

echo "=== [$SP] 4. DIAMOND vs Swiss-Prot ==="
diamond blastp -d "$SPROT" -q ${SP}_cdhit.fasta.transdecoder_dir/longest_orfs.pep \
        --max-target-seqs 1 --evalue 1e-10 --threads "$T" -o blastp.outfmt6

echo "=== [$SP] 5. TransDecoder.Predict ==="
TransDecoder.Predict -t ${SP}_cdhit.fasta \
        --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

echo "=== [$SP] 6. representative isoform per gene ==="
python3 "$SELECT" ${SP}_cdhit.fasta.transdecoder.pep ${SP}_cdhit.fasta.transdecoder.cds \
        ${SP}_longest.isoform.faa ${SP}.cds.fa
echo -n "representative proteins: "; grep -c '^>' ${SP}_longest.isoform.faa
