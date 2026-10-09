#!/usr/bin/env python3
"""Select one representative isoform per Trinity gene.

Rule: the ORF with the highest TransDecoder coding score (score= field in the
TransDecoder.Predict header) per TRINITY_DNx_cx_gx gene; peptide length is used
only as a fallback for headers without a score.  Protein and CDS are written
together so that the two files share sequence identifiers.

The selection is done here in Python on purpose: TransDecoder identifiers contain
'|' characters, which `seqkit grep -r` interprets as a regular-expression
alternation, silently pulling in unrelated sequences.

Usage:
  03_select_representative_isoform.py <pep_in> <cds_in> <pep_out> <cds_out>
"""
import re
import sys

GENE = re.compile(r"(TRINITY_DN\d+_c\d+_g\d+)")
SCORE = re.compile(r"score=(-?[\d.]+)")


def read_fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf)
            name, buf = line[1:].rstrip("\n"), []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    pepsrc, cdssrc, pepout, cdsout = sys.argv[1:5]

    pep = {}
    for h, s in read_fasta(pepsrc):
        pep[h.split()[0]] = (h, s)
    if not pep:
        sys.exit(f"no sequences read from {pepsrc}")

    # rank: coding score when the header carries one, otherwise peptide length
    rank, scored = {}, 0
    for sid, (h, s) in pep.items():
        m = SCORE.search(h)
        if m:
            rank[sid] = float(m.group(1))
            scored += 1
        else:
            rank[sid] = len(s)

    best = {}
    for sid in pep:
        m = GENE.search(sid)
        g = m.group(1) if m else sid
        if g not in best or rank[sid] > rank[best[g]]:
            best[g] = sid
    keep = {best[g] for g in best}

    with open(pepout, "w") as out:
        for g in sorted(best):
            h, s = pep[best[g]]
            out.write(f">{h}\n{s}\n")

    written = 0
    with open(cdsout, "w") as out:
        for h, s in read_fasta(cdssrc):
            sid = h.split()[0]
            if sid in keep:
                out.write(f">{h}\n{s}\n")
                written += 1

    crit = "coding score" if scored == len(pep) else f"mixed ({scored}/{len(pep)} scored)"
    flag = "" if written == len(best) else f"   [WARN cds {written} != pep {len(best)}]"
    print(f"{len(pep):>8} ORFs -> {len(best):>8} genes   (criterion: {crit}){flag}")
    print(f"protein: {pepout}")
    print(f"cds    : {cdsout}")


if __name__ == "__main__":
    main()
