#!/usr/bin/env python3
"""Choose the orthogroup set for the codon analyses and write one CDS file per OG.

Selection: no taxon duplicated (<= 1 copy in each of the seven tree taxa) and at
least MIN_TAXA of the seven present.  Allowing one taxon to be missing keeps every
species in the design while tolerating residual fragmentation of the de novo
transcriptomes; the set stays free of paralogues, so each gene tree remains a
species tree in miniature.  Set A (all seven taxa at exactly one copy; 2,819 OGs)
is recorded in gene_sets.json as a subset flag - this is the set used for the
root-to-tip comparisons (Fig. 1C) and, after QC, for dN/dS (Fig. 1D, E).

Usage:
  02_build_orthogroup_cds.py --ogdir <OrthoFinder Orthogroups dir> \
      --cds <dir with per-species CDS> --out <per-OG dir> --work <dir for gene_sets.json>

Per-species CDS files (label -> file; labels are the tip names used in every tree):
  Colo  Colocasia esculenta      Colo.cds.fa
  La    Landoltia punctata       La.cds.fa
  Laeq  Lemna aequinoctialis     Laeq.cds.fa
  Lgib  Lemna gibba              Lgib.cds.REF_2026-07-22.fa
  Sp    Spirodela polyrhiza      Sp_withorga.cds.renamed.fa
  Wa    Wolffia australiana      Wa_withorga.cds.renamed.fa
  Wo    Wolffiella hyalina       Wh.cds.fa
"""
import os, csv, argparse, json

SPECIES = {                       # OrthoFinder column -> (tip label, cds file)
    "Colocasia_esculenta":  ("Colo", "Colo.cds.fa"),
    "La_longest.isoform":   ("La",   "La.cds.fa"),
    "Laeq_longest.isoform": ("Laeq", "Laeq.cds.fa"),
    "Lgib":                 ("Lgib", "Lgib.cds.REF_2026-07-22.fa"),
    "Sp_withorga":          ("Sp",   "Sp_withorga.cds.renamed.fa"),
    "Wa_withorga":          ("Wa",   "Wa_withorga.cds.renamed.fa"),
    "Wo_longest.isoform":   ("Wo",   "Wh.cds.fa"),
}
MIN_TAXA = 6


def read_fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf)
            name, buf = line[1:].strip().split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ogdir", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    a = ap.parse_args()

    counts = {r["Orthogroup"]: r for r in
              csv.DictReader(open(f"{a.ogdir}/Orthogroups.GeneCount.tsv"), delimiter="\t")}
    genes = {r["Orthogroup"]: r for r in
             csv.DictReader(open(f"{a.ogdir}/Orthogroups.tsv"), delimiter="\t")}
    cols = [c for c in SPECIES if c in next(iter(counts.values()))]
    if len(cols) != 7:
        raise SystemExit(f"expected 7 tree taxa, found {len(cols)}: {cols}")

    def c(og, t):
        return int(counts[og][t])

    allog = list(counts)
    setA = [og for og in allog if all(c(og, t) == 1 for t in cols)]
    print(f"orthogroups total                              {len(allog)}")
    print(f"set A  (7 taxa x exactly 1 copy)               {len(setA)}")
    for k in (7, 6, 5):
        n = sum(1 for og in allog
                if all(c(og, t) <= 1 for t in cols)
                and sum(1 for t in cols if c(og, t) == 1) >= k)
        print(f"set B{k} (no duplicate, >= {k} taxa present)      {n}")
    sel = [og for og in allog
           if all(c(og, t) <= 1 for t in cols)
           and sum(1 for t in cols if c(og, t) == 1) >= MIN_TAXA]
    print(f"\nselected (no duplicate, >= {MIN_TAXA} present):     {len(sel)}")

    seqs = {}
    for col, (lab, f) in SPECIES.items():
        seqs[lab] = dict(read_fasta(f"{a.cds}/{f}"))
        print(f"  loaded {lab:5} {len(seqs[lab]):>7} CDS")

    os.makedirs(a.out, exist_ok=True)
    os.makedirs(a.work, exist_ok=True)
    nw = skipped = 0
    for og in sel:
        row = genes.get(og)
        if row is None:
            skipped += 1
            continue
        recs, ok = [], True
        for col, (lab, _) in SPECIES.items():
            gid = (row.get(col) or "").strip()
            if not gid:
                continue
            gid = gid.split(",")[0].strip()
            s = seqs[lab].get(gid)
            if s is None:
                ok = False
                break
            recs.append((lab, s))
        if not ok or len(recs) < MIN_TAXA:
            skipped += 1
            continue
        os.makedirs(f"{a.out}/{og}", exist_ok=True)
        with open(f"{a.out}/{og}/all.cds.fa", "w") as fh:
            for lab, s in recs:
                fh.write(f">{lab}\n{s}\n")
        nw += 1
    print(f"\nper-OG CDS written: {nw}   skipped: {skipped}")
    json.dump({"selected": sel, "setA": setA}, open(f"{a.work}/gene_sets.json", "w"))
    print(f"gene sets -> {a.work}/gene_sets.json")


if __name__ == "__main__":
    main()
