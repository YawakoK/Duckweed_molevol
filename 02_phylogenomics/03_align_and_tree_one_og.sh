#!/bin/bash
# Codon alignment + maximum-likelihood gene tree for ONE orthogroup.
#
#   1. translate CDS (standard code)
#   2. MAFFT --auto on the proteins
#   3. back-translate the protein alignment to a codon alignment
#   4. IQ-TREE 2 (TIM3+F+I, 1000 ultrafast bootstraps), Colocasia (Colo) as outgroup
#
# Usage:  03_align_and_tree_one_og.sh <OG>        (run in parallel via run_pipeline.sh)
# Env:    TREES   directory holding <OG>/all.cds.fa  (default: TreeEachOG)
#         MAFFT / IQTREE  executables                (default: mafft / iqtree2)
# Versions as run: MAFFT v7.520, IQ-TREE 2.2.6.
set -e
OG=$1
OUT=${TREES:-TreeEachOG}/$OG
MAFFT=${MAFFT:-mafft}
IQ=${IQTREE:-iqtree2}
[ ! -s "$OUT/all.cds.fa" ] && { echo "$OG NOCDS"; exit 0; }
[ -s "$OUT/$OG.treefile" ] && { echo "$OG skip"; exit 0; }

# CDS -> protein (standard code)
python3 - "$OUT/all.cds.fa" "$OUT/all.pep.fa" <<'PY'
import sys
inp,out=sys.argv[1:3]
bases="TCAG"; aas="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codon={}; i=0
for a in bases:
 for b in bases:
  for c in bases:
   codon[a+b+c]=aas[i]; i+=1
def tr(s):
    s=s.upper().replace("U","T"); p=[]
    for k in range(0,len(s)-2,3):
        p.append(codon.get(s[k:k+3],"X"))
    return "".join(p).rstrip("*")
seqs=[];name=None;buf=[]
for L in open(inp):
    if L.startswith(">"):
        if name: seqs.append((name,"".join(buf)))
        name=L[1:].strip();buf=[]
    else: buf.append(L.strip())
if name: seqs.append((name,"".join(buf)))
with open(out,"w") as o:
    for n,s in seqs: o.write(f">{n}\n{tr(s).replace('*','X')}\n")
PY

$MAFFT --auto --quiet "$OUT/all.pep.fa" > "$OUT/all.pep.aln.fa" 2>/dev/null

# back-translate the aligned proteins to a codon alignment
python3 - "$OUT/all.pep.aln.fa" "$OUT/all.cds.fa" "$OUT/all.cds.aln.fa" <<'PY'
import sys
aln,cds,out=sys.argv[1:4]
def rd(f):
    d={};n=None;b=[]
    for L in open(f):
        if L.startswith(">"):
            if n: d[n]="".join(b)
            n=L[1:].strip();b=[]
        else: b.append(L.strip())
    if n: d[n]="".join(b)
    return d
A=rd(aln); C=rd(cds)
with open(out,"w") as o:
    for n,ap in A.items():
        nt=C[n].upper().replace("U","T"); k=0; res=[]
        for aa in ap:
            if aa=="-": res.append("---")
            else:
                res.append(nt[k:k+3] if k+3<=len(nt) else "---"); k+=3
        o.write(f">{n}\n{''.join(res)}\n")
PY

# ML gene tree; if the outgroup is absent from this OG, fall back to an unrooted run
$IQ -s "$OUT/all.cds.aln.fa" -m TIM3+F+I -B 1000 -T 1 -o Colo \
   -pre "$OUT/$OG" -keep-ident -redo -quiet 2>/dev/null || \
$IQ -s "$OUT/all.cds.aln.fa" -m TIM3+F+I -B 1000 -T 1 \
   -pre "$OUT/$OG" -keep-ident -redo -quiet 2>/dev/null
echo "$OG done"
