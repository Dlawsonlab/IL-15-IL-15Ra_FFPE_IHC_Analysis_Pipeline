#!/usr/bin/env python3
"""
Summarise the DAB (IL-15) threshold check run by hpc/run_dab_threshold_check.sh.

Decision rule, fixed before any section result was seen:

  floor(t)       highest DAB-positive area among the four non-DAB controls at threshold t:
                 secondary-only brain, secondary-only liver, haematoxylin-only, Vector Red-only
  separation(t)  each group's median DAB-positive area divided by floor(t)

  Scored on the 21 sections the Vector Red threshold was chosen on.

  0.30 is kept if
    1. every group's median is at least 10x the floor at 0.30, and
    2. the order of the four group medians is the same at every threshold from 0.20 to 0.40.

Also checks reproducibility: the Vector Red scan re-measured in this run must match the
2026-08-12 per-section scans in outputs/threshold_selection/per_sample_scans/.

Usage:  python3 scripts/threshold_evidence/summarize_dab_check.py [--scans DIR]
Writes: outputs/threshold_selection/dab_check/dab_check_summary.csv
        outputs/threshold_selection/dab_check/dab_check_summary.txt
"""
import csv
import glob
import os
import statistics

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
os.chdir(ROOT)
B = "outputs/threshold_selection/dab_check"

# --scans DIR reads the scan CSVs from a copy, for when the mount has not yet shown new files.
# Summaries are always written under B.
import sys
SCANS = sys.argv[sys.argv.index("--scans") + 1] if "--scans" in sys.argv else B
CONTROLS = ["ctrl_secondary_only_brain", "ctrl_secondary_only_liver",
            "ctrl_hematoxylin_only", "ctrl_VectorRed_only"]
GROUP = {"Normal_brain": "Non-neoplastic brain", "BCBM": "BCBM",
         "Normal_liver": "Non-neoplastic liver", "Liver_met": "Liver metastasis"}
ORDER = list(GROUP.values())


def scan(label, dab=True):
    pat = f"{SCANS}/{label}/*_dab_threshold_scan.csv" if dab else f"{SCANS}/{label}/*_threshold_scan.csv"
    files = [f for f in glob.glob(pat) if dab or "_dab_" not in f]
    if not files:
        return None
    col = "dab_area_fraction" if dab else "magenta_area_fraction"
    return {round(float(r["threshold"]), 2): float(r[col]) for r in csv.DictReader(open(files[0]))}


labels = [l.split("|")[0] for l in open("inputs_lists/dab_check_inputs.txt").read().split("\n") if l]
# The same 21 sections the Vector Red threshold was chosen on: those with a 2026-08-12 scan.
# BR-10 is scanned again here but not scored, so both thresholds rest on one section set.
AUG = {os.path.basename(f).split("__")[0] for f in glob.glob("outputs/threshold_selection/per_sample_scans/*_threshold_scan.csv")}
sections = [l for l in labels if not l.startswith("ctrl_") and l in AUG]
ctrl = {c: scan(c) for c in CONTROLS}
sec = {s: scan(s) for s in sections}
missing = [k for k, v in {**ctrl, **sec}.items() if v is None]
done = {s: v for s, v in sec.items() if v is not None}
T = sorted(next(iter(done.values())).keys()) if done else []

rows, lines = [], []
for t in T:
    floor = max(v[t] for v in ctrl.values() if v is not None and t in v)
    med = {}
    for g in ORDER:
        vals = [v[t] for s, v in done.items() if GROUP[s.split("_", 1)[1]] == g]
        med[g] = statistics.median(vals) if vals else float("nan")
    rows.append(dict(threshold=t, floor_pct=floor * 100,
                     **{f"median_pct_{g}": med[g] * 100 for g in ORDER},
                     **{f"separation_{g}": (med[g] / floor if floor > 0 else float("inf")) for g in ORDER}))

os.makedirs(B, exist_ok=True)
if rows:
    with open(f"{B}/dab_check_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

lines.append(f"sections with results: {len(done)}/{len(sections)}; controls: "
             f"{sum(v is not None for v in ctrl.values())}/{len(CONTROLS)}; missing: {missing or 'none'}")
lines.append(f"{'thr':>5} {'floor %':>9} " + " ".join(f"{g[:14]:>15}" for g in ORDER))
for r in rows:
    if r["threshold"] in (0.1, 0.16, 0.2, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.4, 0.5, 0.6):
        lines.append(f"{r['threshold']:5.2f} {r['floor_pct']:9.3f} " +
                     " ".join(f"{r[f'median_pct_{g}']:8.2f}% x{r[f'separation_{g}']:5.0f}" for g in ORDER))

by_t = {r["threshold"]: r for r in rows}
if 0.3 in by_t and not missing:
    r = by_t[0.3]
    c1 = all(r[f"separation_{g}"] >= 10 for g in ORDER)
    rank = lambda t: tuple(sorted(ORDER, key=lambda g: -by_t[t][f"median_pct_{g}"]))
    ref = rank(0.3)
    c2 = all(rank(t) == ref for t in by_t if 0.2 - 1e-9 <= t <= 0.4 + 1e-9)
    lines.append(f"\nrule 1 (every group median >= 10x floor at 0.30): {'PASS' if c1 else 'FAIL'}"
                 f"  -- min separation x{min(r[f'separation_{g}'] for g in ORDER):.0f}")
    lines.append(f"rule 2 (group order unchanged 0.20-0.40): {'PASS' if c2 else 'FAIL'}  -- order {ref}")
    lines.append(f"VERDICT: {'keep 0.30' if c1 and c2 else 'do NOT keep 0.30 without review'}")

# reproducibility: Vector Red re-measured here vs the 2026-08-12 per-section scans
diffs = []
for s in done:
    new = scan(s, dab=False)
    old_f = glob.glob(f"outputs/threshold_selection/per_sample_scans/{s}__*_threshold_scan.csv")
    if new and old_f:
        old = {round(float(x["threshold"]), 2): float(x["magenta_area_fraction"]) for x in csv.DictReader(open(old_f[0]))}
        diffs += [abs(new[t] - old[t]) for t in new if t in old]
if diffs:
    lines.append(f"\nVector Red re-measured vs 2026-08-12 scans: {len(diffs)} values, "
                 f"max |diff| {max(diffs)*100:.5f} percentage points")

open(f"{B}/dab_check_summary.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
