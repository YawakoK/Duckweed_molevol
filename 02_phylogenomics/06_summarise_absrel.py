#!/usr/bin/env python3
"""Per-branch table from the aBSREL JSONs, gene-level QC, and root-to-tip relative
rates from the gene trees.

Writes to --out:
  absrel_per_branch.csv   every branch of every orthogroup (terminal and internal):
                          baseline MG94xREV omega, dS and dN (full adaptive model),
                          aBSREL corrected p, selected-class omega, and passed_QC
  fig1_relative_rates.csv one row per set-A orthogroup (all seven taxa single copy),
                          each column a species' root-to-tip branch length divided
                          by that of Colocasia in the same gene tree  (Fig. 1C)

Gene-level QC (passed_QC = 1), applied to every terminal branch present:
  baseline branch length 0.005-1.0 subs/site; baseline omega <= 2;
  dS 0.01-1.0; no in-frame stop codon and < 80 % gaps in the codon alignment;
  alignment >= 100 codons; and at most two terminal branches with aBSREL
  corrected p < 0.05 (three or more is taken as a sign of alignment error).

Usage:
  06_summarise_absrel.py --json absrel/absrel_json --trees TreeEachOG --out output
"""
import json, glob, os, re, csv, argparse, collections, statistics as st

DUCK = ["Sp", "La", "Laeq", "Lgib", "Wa", "Wo"]
SPP = DUCK + ["Colo", "Pis"]          # Pis (Pistia) is never present in this run
ALPHA = 0.05
STOPS = {"TAA", "TAG", "TGA"}
CLADES = {
    frozenset(DUCK): "Lemnaceae_crown",
    frozenset(["La", "Laeq", "Lgib", "Wa", "Wo"]): "post-Spirodela",
    frozenset(["Laeq", "Lgib", "Wa", "Wo"]): "Lemna+Wolffioideae",
    frozenset(["Laeq", "Lgib"]): "Lemna",
    frozenset(["Wa", "Wo"]): "Wolffioideae",
    frozenset(["La", "Laeq", "Lgib"]): "Landoltia+Lemna",
}


def parse_labelled(s):
    """Newick -> list of (label, tip set, is_leaf, own branch length)."""
    s = s.strip().rstrip(";")
    pos = [0]
    out = []

    def node():
        if s[pos[0]] == "(":
            pos[0] += 1
            kids = []
            while True:
                kids.append(node())
                if s[pos[0]] == ",":
                    pos[0] += 1
                else:
                    break
            pos[0] += 1
            m = re.match(r"[^,;:()]*", s[pos[0]:]); label = m.group(); pos[0] += m.end()
            L = 0.0
            if pos[0] < len(s) and s[pos[0]] == ":":
                pos[0] += 1
                m = re.match(r"[-0-9.eE]+", s[pos[0]:]); L = float(m.group()); pos[0] += m.end()
            tips = set()
            for k in kids:
                tips |= k
            out.append((label, frozenset(tips), False, L))
            return tips
        m = re.match(r"[^,:()]+", s[pos[0]:]); name = m.group(); pos[0] += m.end()
        L = 0.0
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            m = re.match(r"[-0-9.eE]+", s[pos[0]:]); L = float(m.group()); pos[0] += m.end()
        out.append((name, frozenset([name]), True, L))
        return {name}
    node()
    return out


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
    ap.add_argument("--json", required=True, help="directory of aBSREL <OG>.json")
    ap.add_argument("--trees", required=True, help="per-OG directory (all.cds.aln.fa, <OG>.treefile)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    # ------------------------------------------------------------ parse JSONs
    rows = []
    for f in sorted(glob.glob(f"{a.json}/OG*.json")):
        try:
            d = json.load(open(f))
        except Exception:
            continue
        og = os.path.basename(f)[:-5]
        sites = d.get("input", {}).get("number of sites")
        nodes = parse_labelled(d["input"]["trees"]["0"])
        tipmap = {l: (ts, leaf) for l, ts, leaf, _ in nodes}
        for b, at in d["branch attributes"]["0"].items():
            ts, leaf = tipmap.get(b, (frozenset(), True))
            if leaf:
                label = b if b in SPP else None
            else:
                duck = frozenset(ts & set(DUCK))
                label = CLADES.get(duck) if duck and not (ts & {"Colo", "Pis"}) else None
            rd = at.get("Rate Distributions") or []
            pos = [(w, fr) for w, fr in rd if w > 1]
            rows.append(dict(OG=og, branch=b, leaf=leaf, label=label, sites=sites,
                             bl=at.get("Baseline MG94xREV"),
                             om=at.get("Baseline MG94xREV omega ratio"),
                             dS=at.get("Full adaptive model (synonymous subs/site)"),
                             dN=at.get("Full adaptive model (non-synonymous subs/site)"),
                             p=at.get("Corrected P-value"), LRT=at.get("LRT"),
                             w_sel=max((w for w, _ in pos), default=None),
                             f_sel=sum(fr for _, fr in pos) if pos else 0.0))
    by_og = collections.defaultdict(dict)
    for r in rows:
        by_og[r["OG"]][r["branch"]] = r
    ALL = sorted(by_og)
    print(f"orthogroups with aBSREL output : {len(ALL)}")
    print(f"branches parsed                : {len(rows)} "
          f"({sum(1 for r in rows if r['leaf'])} terminal / {sum(1 for r in rows if not r['leaf'])} internal)")

    # ------------------------------------------------------------ QC
    stopinfo = collections.defaultdict(dict)
    for og in ALL:
        f = f"{a.trees}/{og}/all.cds.aln.fa"
        if not os.path.exists(f):
            continue
        for sp, s in read_fasta(f):
            if sp not in SPP:
                continue
            s = s.upper().replace("U", "T")
            L = (len(s) // 3) * 3
            stopinfo[og][sp] = (sum(1 for i in range(0, L, 3) if s[i:i+3] in STOPS),
                                s.count("-") / max(len(s), 1))

    def gene_ok(og):
        d = by_og[og]
        present = [sp for sp in SPP if sp in d]
        if len(present) < 6:
            return False
        for sp in present:
            r = d[sp]
            if r["bl"] is None or not (0.005 <= r["bl"] <= 1.0):
                return False
            if r["om"] is None or r["om"] > 2.0:
                return False
            if r["dS"] is None or not (0.01 <= r["dS"] <= 1.0):
                return False
            si = stopinfo.get(og, {}).get(sp)
            if si is None or si[0] > 0 or si[1] > 0.8:
                return False
        ref = next(iter(d.values()))
        return bool(ref["sites"] and ref["sites"] >= 100)

    KEEP = [og for og in ALL if gene_ok(og)]
    multi = {og for og in KEEP
             if sum(1 for sp in SPP if by_og[og].get(sp) and by_og[og][sp]["p"] is not None
                    and by_og[og][sp]["p"] < ALPHA) >= 3}
    CLEAN = set(og for og in KEEP if og not in multi)
    print(f"gene-level QC passed           : {len(KEEP)}")
    print(f"after dropping multi-flagged   : {len(CLEAN)}  (dropped {len(multi)})")

    with open(f"{a.out}/absrel_per_branch.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["OG", "branch", "is_terminal", "lineage", "baseline_omega", "corrected_p",
                    "LRT", "omega_selected", "prop_sites_selected", "dS", "dN", "sites",
                    "passed_QC"])
        for og in ALL:
            for b, r in by_og[og].items():
                w.writerow([og, b, int(r["leaf"]), r["label"] or "", r["om"], r["p"], r["LRT"],
                            r["w_sel"], round(r["f_sel"], 4), r["dS"], r["dN"], r["sites"],
                            int(og in CLEAN)])

    # ------------------------------------------------------------ Fig. 1C quantity
    # Root-to-tip branch length per duckweed, divided by that of Colocasia in the
    # same gene tree.  The root is the Colocasia / Lemnoideae split; a tip's path
    # sums its own edge and every internal edge whose tip set contains it (the
    # root itself, which subtends all seven taxa, is excluded).
    print("\nroot-to-tip relative rate per duckweed (normalised by Colocasia)")
    rate = collections.defaultdict(list)
    fig = []
    for d in sorted(glob.glob(f"{a.trees}/OG*")):
        og = os.path.basename(d)
        tf = f"{d}/{og}.treefile"
        if not os.path.exists(tf):
            continue
        try:
            nodes = parse_labelled(open(tf).read())
        except Exception:
            continue
        lens = {lab: L for lab, ts, leaf, L in nodes}
        inner = [(lab, ts, L) for lab, ts, leaf, L in nodes if not leaf]

        def root_to_tip(sp):
            tot = lens.get(sp, 0.0)
            for lab, ts, L in inner:
                if sp in ts and len(ts) < 7:
                    tot += L
            return tot
        if "Colo" not in lens or root_to_tip("Colo") <= 0:
            continue
        base = root_to_tip("Colo")
        row = {"OG": og}
        ok = True
        for sp in DUCK:
            if sp not in lens:
                ok = False
                break
            row[sp] = root_to_tip(sp) / base
        if not ok:
            continue
        fig.append(row)
        for sp in DUCK:
            rate[sp].append(row[sp])
    print(f"  gene trees used: {len(fig)}")
    for sp in sorted(DUCK, key=lambda x: st.median(rate[x])):
        print(f"    {sp:6} median {st.median(rate[sp]):.3f}")
    with open(f"{a.out}/fig1_relative_rates.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["OG"] + DUCK)
        for row in fig:
            w.writerow([row["OG"]] + [f"{row[sp]:.6f}" for sp in DUCK])
    print(f"\nwrote {a.out}/absrel_per_branch.csv and {a.out}/fig1_relative_rates.csv")


if __name__ == "__main__":
    main()
