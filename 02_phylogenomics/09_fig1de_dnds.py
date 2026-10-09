#!/usr/bin/env python3
"""Fig. 1D, 1E: synonymous and non-synonymous divergence per lineage.

Inputs (from data/ by default; set IN=output to use regenerated tables):
  absrel_per_branch.csv   per-branch dN, dS, selected-class omega, passed_QC (06)
  roottotip_dnds.csv      root-to-tip dS/dN per species per orthogroup       (07)
  Orthogroups.tsv         OrthoFinder table, to restrict to single-copy OGs

Orthogroups used = passed_QC AND single copy in all seven non-Arabidopsis taxa
(2,474; main text L176).

Panel D: root-to-tip dS per species (violins + boxes).
Panel E: median terminal-branch dN against median dS per species (bars = IQR);
         terminal branches with a selected-class omega > 5 are excluded.  The grey
         line is the proportional fit through the origin (constant omega); its
         slope and the per-species omega are printed (main text L181-183).

Output: output/Fig1DE_dN_dS.pdf/.png
"""
import os, csv, collections, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

IN = os.environ.get("IN", "data")
OUT = "output"
os.makedirs(OUT, exist_ok=True)
PB = f"{IN}/absrel_per_branch.csv"
RTT = f"{IN}/roottotip_dnds.csv"
OGTSV = "data/Orthogroups.tsv"
OUTFILE = f"{OUT}/Fig1DE_dN_dS"

ANTHO = {"Sp": True, "La": True, "Lgib": True,
         "Laeq": False, "Wa": False, "Wo": False}

# per-species colour and marker, matching the shared species legend
COL = {"Sp": "#4a3d70", "La": "#6b5c96", "Lgib": "#a99cc9",
       "Laeq": "#8fb85a", "Wo": "#3f5d24", "Wa": "#6f9a3c"}
MRK = {"Sp": "o", "La": "^", "Lgib": "s",
       "Laeq": "D", "Wo": (4, 2, 0), "Wa": (6, 2, 0)}   # plus / asterisk are stroked
STROKE = {"Wo", "Wa"}
ABBR = {"Sp": "Sp", "La": "Lp", "Lgib": "Lg",
        "Laeq": "La", "Wa": "Wa", "Wo": "Wh"}
FULLNAME = {"Sp": "Spirodela polyrhiza", "La": "Landoltia punctata",
            "Lgib": "Lemna gibba", "Laeq": "Lemna aequinoctialis",
            "Wo": "Wolffiella hyalina", "Wa": "Wolffia australiana"}


def qc_ogs():
    b = collections.defaultdict(list)
    for r in csv.DictReader(open(PB)):
        b[r["OG"]].append(r)
    return {og for og, rs in b.items() if all(r["passed_QC"] == "1" for r in rs)}


def single_copy_ogs():
    out = set()
    with open(OGTSV) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        idx = [i for i, h in enumerate(hdr)
               if h not in ("Orthogroup", "Araport11_pep_20250214")]
        for row in rd:
            if all(row[i].strip() and "," not in row[i] for i in idx):
                out.add(row[0])
    return out


KEEP = qc_ogs() & single_copy_ogs()
print("orthogroups used:", len(KEEP))          # -> 2474 (main text L176)

rtt = collections.defaultdict(list)
for r in csv.DictReader(open(RTT)):
    if r["OG"] in KEEP:
        try:
            rtt[r["species"]].append(float(r["dS_total"]))
        except ValueError:
            pass

tdn, tds = collections.defaultdict(list), collections.defaultdict(list)
for r in csv.DictReader(open(PB)):
    if r["OG"] not in KEEP or r["is_terminal"] not in ("1", "True", "TRUE"):
        continue
    lin = r["lineage"]
    if lin not in ANTHO:
        continue
    ws = r.get("omega_selected", "")
    if ws not in ("", "NA", "None"):
        try:
            if float(ws) > 5:
                continue
        except ValueError:
            pass
    try:
        tdn[lin].append(float(r["dN"]))
        tds[lin].append(float(r["dS"]))
    except ValueError:
        pass

order = [s for s in ["Sp", "La", "Lgib", "Laeq", "Wa", "Wo"] if s in rtt]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(6.5, 3.1))

# ---------------- Panel D ----------------
data = [rtt[s] for s in order]
vp = axA.violinplot(data, positions=range(len(order)), widths=0.8,
                    showextrema=False, showmedians=False)
for body, s in zip(vp["bodies"], order):
    body.set_facecolor(COL[s]); body.set_edgecolor("gray")
    body.set_alpha(1.0); body.set_linewidth(0.4)
axA.boxplot(data, positions=range(len(order)), widths=0.13,
            showfliers=False, patch_artist=True,
            boxprops=dict(facecolor="white", linewidth=0.6),
            medianprops=dict(color="black", linewidth=0.8),
            whiskerprops=dict(linewidth=0.6), capprops=dict(linewidth=0.6))
axA.set_xticks(range(len(order)))
axA.set_xticklabels([FULLNAME[s] for s in order], fontsize=6.5,
                    style="italic", rotation=40, ha="right",
                    rotation_mode="anchor")
axA.set_ylabel("Root-to-tip dS", fontsize=8)
axA.tick_params(labelsize=7)
axA.set_ylim(0, 0.55)
axA.grid(axis="y", color="0.9", linewidth=0.5); axA.set_axisbelow(True)
handles = [plt.Rectangle((0, 0), 1, 1, fc="#7a6aa8", ec="gray", lw=0.4),
           plt.Rectangle((0, 0), 1, 1, fc="#6f9a3c", ec="gray", lw=0.4)]
axA.legend(handles, ["Antho. accum.", "Antho. non-accum."],
           fontsize=6.5, loc="upper left", frameon=True)
axA.text(-0.16, 1.04, "D", transform=axA.transAxes,
         fontsize=12, fontweight="bold", va="top")
print("median root-to-tip dS:", {s: round(st.median(rtt[s]), 3) for s in order})

# ---------------- Panel E ----------------
xs, ys = [], []
sp_handles = []
for s in order:
    mds, mdn = st.median(tds[s]), st.median(tdn[s])
    xs.append(mds); ys.append(mdn)
    qs = np.percentile(tds[s], [25, 75]); qn = np.percentile(tdn[s], [25, 75])
    if s in STROKE:                       # plus / asterisk: drawn as strokes
        mstyle = dict(mfc="none", mec=COL[s], mew=2.0, ms=9.5)
    else:
        mstyle = dict(mfc=COL[s], mec="black", mew=0.6, ms=6)
    axB.errorbar(mds, mdn,
                 xerr=[[mds - qs[0]], [qs[1] - mds]],
                 yerr=[[mdn - qn[0]], [qn[1] - mdn]],
                 fmt="none", ecolor=COL[s], elinewidth=0.6, capsize=1.5,
                 alpha=0.4, zorder=2)
    h, = axB.plot([mds], [mdn], marker=MRK[s], linestyle="none",
                  label=f"$\\it{{{ABBR[s]}}}$  $\\omega$={mdn/mds:.2f}",
                  zorder=4, **mstyle)
    sp_handles.append(h)
    print(f"  {ABBR[s]:3} median dS {mds:.4f}  median dN {mdn:.4f}  omega {mdn/mds:.2f}")

omega = sum(x * y for x, y in zip(xs, ys)) / sum(x * x for x in xs)
xlim = max(xs) * 1.28
line, = axB.plot([0, xlim], [0, omega * xlim], ls="-", color="0.68", lw=0.9,
                 label=f"proportional ($\\omega$={omega:.2f})", zorder=1)
axB.set_xlim(0, xlim); axB.set_ylim(0, max(ys) * 1.9)
axB.set_xlabel("Median dS per lineage ($\\pm$IQR)", fontsize=8)
axB.set_ylabel("Median dN per lineage ($\\pm$IQR)", fontsize=8)
axB.tick_params(labelsize=7)
axB.legend(handles=[line] + sp_handles, fontsize=5.5, loc="upper left",
           ncol=2, frameon=True, handletextpad=0.4, labelspacing=0.3,
           columnspacing=0.9, borderpad=0.4, borderaxespad=0.4)
axB.grid(color="0.9", linewidth=0.5); axB.set_axisbelow(True)
axB.text(-0.18, 1.04, "E", transform=axB.transAxes,
         fontsize=12, fontweight="bold", va="top")

print(f"proportional fit through origin: omega = {omega:.3f}")   # -> 0.50 (main text L181)
fig.tight_layout()
fig.savefig(OUTFILE + ".pdf")
fig.savefig(OUTFILE + ".png", dpi=600)
print("wrote", OUTFILE + ".pdf/.png")
