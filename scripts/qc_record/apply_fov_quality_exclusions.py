#!/usr/bin/env python3
"""
Apply the reviewer-approved quality exclusions to the per-sample *_fov_exclude.csv files.

The reviewer approved the full flagged set in outputs/tables/fov_quality_flagged.csv on
2026-08-12, working from outputs/reports/fov_quality_review.html. ANALYSIS_LOG §29 lists the
resulting exclusions by verdict. The script writes each flagged field into its sample's
exclusion file with its verdict and the metrics behind it. Every exclusion therefore traces
to a number and a rule.

Existing exclusions are preserved and never rewritten, so earlier reviewer decisions keep
their recorded reasons. The script skips any field that is already excluded, so re-running
it is safe.

Run with --dry-run to see what would change without writing.
"""
import csv, glob, os, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
DRY = "--dry-run" in sys.argv

coh = {r["sample_id"]: r for r in csv.DictReader(open("outputs/tables/cohort_table.csv"))}
flagged = list(csv.DictReader(open("outputs/tables/fov_quality_flagged.csv")))

by_sample = defaultdict(list)
for r in flagged:
    by_sample[r["sample_id"]].append(r)

added_total, skipped_total, files = 0, 0, 0
for sid, rows in sorted(by_sample.items()):
    sub = coh.get(sid, {}).get("subfolder")
    if not sub:
        print(f"  {sid}: not in cohort table, skipped"); continue
    path = os.path.join("per_sample_out", sub, f"{sub}_fov_exclude.csv")

    existing = {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            existing[str(int(r["fov_id"]))] = r.get("reason", "")

    added = 0
    for r in rows:
        fid = str(int(r["fov_id"]))
        if fid in existing:
            skipped_total += 1
            continue
        reason = (f"reviewer-approved 2026-08-12: quality flag [{r['verdict']}] at native "
                  f"resolution (void {r['void_total']}, largest {r['tear_area']}, "
                  f"solidity {r['tear_solidity']}, dark {r['fold_area']}, "
                  f"elong {r['fold_elongation']})")
        existing[fid] = reason
        added += 1

    if added:
        files += 1
        added_total += added
        if not DRY:
            with open(path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["fov_id", "reason"]); w.writeheader()
                for fid in sorted(existing, key=int):
                    w.writerow({"fov_id": fid, "reason": existing[fid]})
        print(f"  {sid:13s} +{added:3d} excluded  (file now holds {len(existing)})")

print(f"\n  {'DRY RUN -- ' if DRY else ''}{added_total} field(s) newly excluded across {files} sample(s); "
      f"{skipped_total} already excluded and left untouched")
