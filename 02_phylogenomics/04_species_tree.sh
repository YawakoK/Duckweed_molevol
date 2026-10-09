#!/bin/bash
# Concatenated species tree (Fig. 1B).
#
#   single-copy codon alignments (7 taxa, one sequence each)
#     -> AMAS concat (one partition per orthogroup)
#     -> IQ-TREE 3 partitioned ML, 1000 ultrafast bootstraps, Colocasia outgroup
#
# The tree topology and support values were then drawn with PhyloWeaver
# (https://yawak.jp/PhyloWeaver/) for the figure.
#
# Env:  TREES   per-OG directory from 03_align_and_tree_one_og.sh (default: TreeEachOG)
#       OUT     output directory (default: species_tree)
#       AMAS    path to AMAS.py ; IQTREE3 executable (default: iqtree3) ; NT threads
# Versions as run: AMAS (Borowiec 2016), IQ-TREE 3.0.1.
set -e
TREES=${TREES:-TreeEachOG}
OUT=${OUT:-species_tree}
AMAS=${AMAS:-AMAS.py}
IQ=${IQTREE3:-iqtree3}
NT=${NT:-24}

mkdir -p "$OUT"
cd "$OUT"

# 1. collect single-copy codon alignments (exactly 7 distinct taxa)
: > sc_alignments.txt
for d in "$TREES"/OG*; do
    a="$d/all.cds.aln.fa"
    [ -s "$a" ] || continue
    n=$(grep -c '^>' "$a")
    [ "$n" -eq 7 ] || continue
    u=$(grep '^>' "$a" | sed 's/^>//' | sort -u | wc -l)
    [ "$u" -eq 7 ] || continue
    echo "$a" >> sc_alignments.txt
done
echo "single-copy alignments: $(wc -l < sc_alignments.txt)"

# 2. concatenate, one partition per orthogroup
python3 "$AMAS" concat -f fasta -d dna \
        -i $(cat sc_alignments.txt) \
        -p partition.txt -t Concatenated.fna
awk '{print "DNA, " $1 " = " $3}' partition.txt > partition_iqtree.txt

# 3. partitioned ML tree
$IQ -s Concatenated.fna -p partition_iqtree.txt \
    -B 1000 -T "$NT" -redo -o Colo -pre species_tree > iqtree.out 2>&1

cat species_tree.treefile
