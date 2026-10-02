#!/usr/bin/env python3
"""
Summarise the Vector Red (IL-15Ra) threshold scan run by hpc/run_il15_threshold_scan.sh
(job 55210127, 2026-08-12), the evidence behind MAGENTA_THRESHOLD = 0.10.

  evidence   every section x threshold: positive-area fraction and area, stacked from the
             per-section scan CSVs
  summary    per threshold, each group's median positive area (% of tissue, 4 decimals) and
             the liver-to-brain ratio of those rounded medians

The scan covered the 22 sections in inputs_lists/thrscan_samples.txt. One task wrote a 0-byte
CSV (BR-10, normal brain), so 21 are summarised. The script checks the count against the
input list and names any section that is missing, rather than summing whatever files exist.

Usage:  python3 scripts/threshold_evidence/summarize_threshold_scan.py
Writes: outputs/threshold_selection/threshold_scan_evidence.csv
        outputs/threshold_selection/threshold_scan_summary.csv
"""
import csv
import glob
import os
import statistics

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
os.chdir(ROOT)
B = "outputs/threshold_selection"
GROUP = {"Normal_brain": "brain", "BCBM": "met", "Normal_liver": "liver", "Liver_met": "liver_met"}
KNOWN_MISSING = {"BR-10"}       # its scan task wrote a 0-byte CSV

subdirs = open("subdirs.txt").read().splitlines()
expected = {subdirs[int(i) - 1].split("_")[0] for i in open("inputs_lists/thrscan_samples.txt").read().split()}

rows = []
for p in sorted(glob.glob(f"{B}/per_sample_scans/*_threshold_scan.csv")):
    sub = os.path.basename(p).split("__")[0]
    sid = sub.split("_")[0]
    grp = next(g for k, g in GROUP.items() if sub.endswith(k))
    for r in csv.DictReader(open(p)):
        rows.append({"sample_id": sid, "group": grp, "threshold": round(float(r["threshold"]), 2),
                     "magenta_area_fraction": float(r["magenta_area_fraction"]),
                     "positive_area_mm2": float(r["positive_area_mm2"])})

found = {r["sample_id"] for r in rows}
missing = expected - found
assert found <= expected, f"scan CSVs for sections not in the input list: {sorted(found - expected)}"
assert missing == KNOWN_MISSING, f"missing scans differ from the known gap: {sorted(missing)}"
print(f"  {len(found)} of {len(expected)} sections scanned; missing {sorted(missing)} (known)")

with open(f"{B}/threshold_scan_evidence.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

thresholds = sorted({r["threshold"] for r in rows})
with open(f"{B}/threshold_scan_summary.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["threshold", "brain_median_pct", "met_median_pct", "liver_median_pct",
                "liver_met_median_pct", "liver_to_brain_ratio"])
    for t in thresholds:
        med = {g: round(statistics.median(r["magenta_area_fraction"] * 100 for r in rows
                                          if r["threshold"] == t and r["group"] == g), 4)
               for g in GROUP.values()}
        w.writerow([f"{t:.2f}", f"{med['brain']:.5f}", f"{med['met']:.4f}", f"{med['liver']:.4f}",
                    f"{med['liver_met']:.4f}", round(med["liver"] / med["brain"])])
print(f"  wrote {B}/threshold_scan_evidence.csv ({len(rows)} rows), "
      f"{B}/threshold_scan_summary.csv ({len(thresholds)} thresholds)")
