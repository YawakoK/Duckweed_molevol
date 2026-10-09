#!/usr/bin/env python3
"""Root-to-tip dS and dN accumulated from the duckweed common ancestor (Fig. 1D).

For every QC-passing orthogroup, per-branch dN and dS (aBSREL full adaptive model,
as parsed into absrel_per_branch.csv) are summed along the path from the duckweed
MRCA to each duckweed tip; the Colocasia edge is excluded.  Because dN on branches
with an omega > 1 rate class is somewhat inflated, root-to-tip omega on lineages
with more positive selection is a slight over-estimate; only dS is plotted.

Usage:
  07_roottotip_dnds.py --json absrel/absrel_json --per-branch output/absrel_per_branch.csv --out output
Writes: <out>/roottotip_dnds.csv  (OG, species, dN_total, dS_total, omega_rootToTip)
"""
import json, glob, os, re, csv, argparse, collections, statistics as st

DUCK = ["Sp", "La", "Laeq", "Lgib", "Wa", "Wo"]


def parse(tree):
    s = tree.strip().rstrip(";")
    pos = [0]
    nodes = []

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
            m = re.match(r"[^,;:()]*", s[pos[0]:]); name = m.group(); pos[0] += m.end()
            if pos[0] < len(s) and s[pos[0]] == ":":
                pos[0] += 1
                m = re.match(r"[-0-9.eE]+", s[pos[0]:]); pos[0] += m.end()
            idx = len(nodes)
            nodes.append(dict(name=name, children=kids, leaf=False))
            return idx
        m = re.match(r"[^,:()]+", s[pos[0]:]); name = m.group(); pos[0] += m.end()
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            m = re.match(r"[-0-9.eE]+", s[pos[0]:]); pos[0] += m.end()
        idx = len(nodes)
        nodes.append(dict(name=name, children=[], leaf=True))
        return idx
    root = node()
    return nodes, root


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--per-branch", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    QC = {r["OG"] for r in csv.DictReader(open(a.per_branch)) if r["passed_QC"] == "1"}

    out_rows = []
    n_used = 0
    for f in sorted(glob.glob(f"{a.json}/OG*.json")):
        og = os.path.basename(f)[:-5]
        if og not in QC:
            continue
        try:
            d = json.load(open(f))
        except Exception:
            continue
        ba = d["branch attributes"]["0"]

        def dn(name):
            return ba.get(name, {}).get("Full adaptive model (non-synonymous subs/site)")

        def ds(name):
            return ba.get(name, {}).get("Full adaptive model (synonymous subs/site)")
        nodes, root = parse(d["input"]["trees"]["0"])

        def walk(idx, accN, accS):
            nd = nodes[idx]
            nm = nd["name"]
            if idx != root:
                x = dn(nm); y = ds(nm)
                if x is None or y is None:
                    return
                accN += x; accS += y
            if nd["leaf"]:
                if nm in DUCK and accS and accS > 0:
                    out_rows.append((og, nm, round(accN, 5), round(accS, 5), round(accN / accS, 5)))
                return
            for c in nd["children"]:
                if nodes[c]["leaf"] and nodes[c]["name"] == "Colo":
                    continue
                walk(c, accN, accS)
        walk(root, 0.0, 0.0)
        n_used += 1

    with open(f"{a.out}/roottotip_dnds.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["OG", "species", "dN_total", "dS_total", "omega_rootToTip"])
        w.writerows(out_rows)
    print(f"gene trees used: {n_used}")
    print(f"rows (species x gene): {len(out_rows)} -> {a.out}/roottotip_dnds.csv")
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for og, sp, x, y, w_ in out_rows:
        by[sp]["dN"].append(x); by[sp]["dS"].append(y); by[sp]["om"].append(w_)
    print(f"\n{'sp':6}{'n':>6}{'dS med':>9}{'dN med':>9}{'w med':>9}")
    for sp in DUCK:
        v = by[sp]
        print(f"{sp:6}{len(v['dN']):>6}{st.median(v['dS']):>9.3f}{st.median(v['dN']):>9.3f}{st.median(v['om']):>9.3f}")


if __name__ == "__main__":
    main()
