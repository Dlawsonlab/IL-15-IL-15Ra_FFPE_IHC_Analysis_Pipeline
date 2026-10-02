#!/usr/bin/env python3
"""
Test whether the IL-15 (DAB) group contrasts depend on the positivity threshold.

Endpoint: whole-slide DAB-positive area of all 57 analysed sections, from the threshold scans
of hpc/run_dab_threshold_check.sh. The 19 default-geometry sections come from the check run and the
other 38 from the endpoint run. Each section's value at 0.30 must reproduce whole-slide IL-15
in cohort_table.csv, or the script stops.

Statistics follow METHODS section 10: Kruskal-Wallis; two-sided Mann-Whitney U on all six
pairs; Benjamini-Hochberg across the six.

Criterion (set before results): every contrast significant at 0.30 stays
significant with the same direction at every threshold from 0.24 to 0.40, and no contrast
non-significant at 0.30 becomes significant.

Usage:  python3 scripts/threshold_evidence/dab_threshold_sensitivity.py [--check DIR] [--endpoint DIR]
        (DIRs hold the scan folders; they default to outputs/threshold_selection/dab_check
        and dab_endpoint, and can point at copies when the mount lags)
Writes: outputs/threshold_selection/dab_endpoint/dab_sensitivity.csv
        outputs/threshold_selection/dab_endpoint/dab_sensitivity.txt
"""
import csv
import glob
import os
import sys

import numpy as np
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
os.chdir(ROOT)
T = "outputs/threshold_selection"


def arg(flag, default):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else default


CHECK, ENDPOINT = arg("--check", f"{T}/dab_check"), arg("--endpoint", f"{T}/dab_endpoint")
GROUPS = ["Normal brain", "BCBM", "Normal liver", "Liver met"]
LABEL = {"Normal brain": "non-neoplastic brain", "BCBM": "BCBM",
         "Normal liver": "non-neoplastic liver", "Liver met": "liver metastasis"}
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
WINDOW = [round(0.24 + 0.02 * k, 2) for k in range(9)]          # 0.24 ... 0.40


def bh(p):
    m = len(p); order = np.argsort(p); adj = [0.0] * m; prev = 1.0
    for k, i in enumerate(reversed(order)):
        prev = min(prev, p[i] * m / (m - k)); adj[i] = prev
    return adj


def scan(folder):
    f = glob.glob(f"{folder}/*_dab_threshold_scan.csv")
    return {round(float(r["threshold"]), 2): float(r["dab_area_fraction"])
            for r in csv.DictReader(open(f[0]))} if f else None


coh = list(csv.DictReader(open("outputs/tables/cohort_table.csv")))
endpoint_set = {l.split("|")[0] for l in open("inputs_lists/dab_endpoint_inputs.txt").read().split("\n") if l}
vals, missing, worst = {}, [], 0.0
for r in coh:
    sub = r["subfolder"]
    s = scan(f"{ENDPOINT}/{sub}" if sub in endpoint_set else f"{CHECK}/{sub}")
    if s is None:
        missing.append(sub); continue
    worst = max(worst, abs(s[0.3] - float(r["ws_brown_af"])))
    vals[sub] = (r["tissue_type"], s)
if missing:
    sys.exit(f"missing scans for {len(missing)} section(s): {missing}")
if worst > 1e-6:
    sys.exit(f"scan does not reproduce whole-slide IL-15 at 0.30 (max diff {worst:.2e}); stopping")

rows, out = [], []
out.append(f"sections: {len(vals)}; scan reproduces whole-slide IL-15 at 0.30 "
           f"(max |diff| {worst * 100:.6f} percentage points)")
for t in sorted({*WINDOW, 0.30}):
    g = [np.array([s[t] for tt, s in vals.values() if tt == grp]) * 100 for grp in GROUPS]
    kw = stats.kruskal(*g).pvalue
    raw = [stats.mannwhitneyu(g[i], g[j], alternative="two-sided").pvalue for i, j in PAIRS]
    adj = bh(raw)
    for (i, j), r_, a in zip(PAIRS, raw, adj):
        rows.append(dict(threshold=t, pair=f"{LABEL[GROUPS[i]]} vs {LABEL[GROUPS[j]]}",
                         median_a=float(np.median(g[i])), median_b=float(np.median(g[j])),
                         direction=">" if np.median(g[i]) > np.median(g[j]) else "<",
                         p=r_, p_adj=a, kruskal_p=kw))

os.makedirs(f"{T}/dab_endpoint", exist_ok=True)
with open(f"{T}/dab_endpoint/dab_sensitivity.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

pairs = [f"{LABEL[GROUPS[i]]} vs {LABEL[GROUPS[j]]}" for i, j in PAIRS]
out.append("\nBH-adjusted p by threshold (direction: > or < between the pair's medians)")
out.append(f"{'pair':44s} " + " ".join(f"{t:>9.2f}" for t in WINDOW))
for pr in pairs:
    cells = []
    for t in WINDOW:
        r = next(x for x in rows if x["threshold"] == t and x["pair"] == pr)
        cells.append(f"{r['direction']}{r['p_adj']:.4f}{'*' if r['p_adj'] < 0.05 else ' '}")
    out.append(f"{pr:44s} " + " ".join(f"{c:>9s}" for c in cells))

ref = {r["pair"]: r for r in rows if r["threshold"] == 0.30}
fails = []
for t in WINDOW:
    for pr in pairs:
        r, r0 = next(x for x in rows if x["threshold"] == t and x["pair"] == pr), ref[pr]
        if r0["p_adj"] < 0.05 and (r["p_adj"] >= 0.05 or r["direction"] != r0["direction"]):
            fails.append(f"{pr}: significant at 0.30, not at {t:.2f}")
        if r0["p_adj"] >= 0.05 and r["p_adj"] < 0.05:
            fails.append(f"{pr}: non-significant at 0.30, significant at {t:.2f}")
out.append("\nVERDICT: " + ("ROBUST -- every contrast keeps its status and direction from 0.24 to 0.40"
                            if not fails else "NOT ROBUST -- " + "; ".join(fails)))
open(f"{T}/dab_endpoint/dab_sensitivity.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
