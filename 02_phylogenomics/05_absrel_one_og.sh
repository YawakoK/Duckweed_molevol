#!/bin/bash
# HyPhy aBSREL on ONE orthogroup (all branches tested).
#
# aBSREL first fits the MG94xREV codon model to every branch (baseline branch
# length and omega) and then the adaptive branch-site model; both sets of
# per-branch estimates are read from the JSON by 06_summarise_absrel.py.
# In-frame stop codons are masked (---) before the fit.
#
# Usage:  05_absrel_one_og.sh <OG>
# Env:    TREES  per-OG directory (default: TreeEachOG)
#         ABSREL output directory for JSONs (default: absrel/absrel_json)
#         HYPHY  executable (default: hyphy)
# Version as run: HyPhy 2.5.61.
set -e
OG=$1
SRC=${TREES:-TreeEachOG}/$OG
OUTDIR=${ABSREL:-absrel/absrel_json}
OUT=$OUTDIR/$OG.json
HY=${HYPHY:-hyphy}
mkdir -p "$OUTDIR"
[ -s "$OUT" ] && exit 0
[ ! -f "$SRC/all.cds.aln.fa" ] && exit 0
[ ! -f "$SRC/$OG.treefile"   ] && exit 0
WORK=$(mktemp -d -t absrel_${OG}_XXXX)
python3 - "$SRC/all.cds.aln.fa" "$WORK/aln.fa" <<'PY'
import sys
inp,out=sys.argv[1:3]
stops={"TAA","TAG","TGA"}
seqs={}; n=None; b=[]
for L in open(inp):
    if L.startswith(">"):
        if n: seqs[n]="".join(b)
        n=L[1:].strip(); b=[]
    else: b.append(L.strip())
if n: seqs[n]="".join(b)
with open(out,"w") as o:
    for n,s in seqs.items():
        s=s.upper().replace("U","T")
        out_s=[]
        L=(len(s)//3)*3
        for i in range(0,L,3):
            c=s[i:i+3]
            out_s.append("---" if c in stops else c)
        o.write(f">{n}\n{''.join(out_s)}\n")
PY
cp "$SRC/$OG.treefile" "$WORK/tree.nwk"
OUT_ABS=$(readlink -f "$OUT")
cd "$WORK"
$HY CPU=1 absrel --alignment aln.fa --tree tree.nwk --output "$OUT_ABS" --branches All > absrel.log 2>&1 || {
  echo "$OG FAIL"; head -3 errors.log 2>/dev/null; }
cd /; rm -rf "$WORK"
echo "$OG done"
