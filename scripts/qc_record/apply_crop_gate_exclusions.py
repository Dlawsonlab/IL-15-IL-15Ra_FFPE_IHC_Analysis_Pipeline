#!/usr/bin/env python3
"""
Apply the reviewer-approved CROP GATE to the per-sample exclusion lists.

The reviewer set these rules on the round-2 review pages. ANALYSIS_LOG §31 records the
evidence behind them; outputs/reports/crop_gate_review.html lists the fields kept by adjudication.

  1. VOID rule. A field is damaged when its largest single near-white component covers
     >= 2% of the field. Such fields are ragged, torn, or straddle the section margin.
     Ragged fields score 0.047-0.111 and solid fields 0.001-0.006, so the cut has an ~8x
     margin. The rule measures the largest void because a total-void rule would discard
     steatosis and liver sinusoids and so bias the liver groups.

  2. SPARSE rule. A field fails when its tissue fraction is <= 80%.

  3. DARK rule. A field fails when its dark fraction is >= 1% of the field. The reviewer
     declined a connectivity rule that separates folds (one large elongated dark region)
     from dense chromatin (thousands of tiny ones). Every field in the dark-only review is
     excluded, apart from the LM-09 fields in rule 4. The connectivity metrics stay in
     outputs/tables/crop_scores_all.json for reference, and they gate nothing.

  4. LM-09 is adjudicated field by field and skips the automatic rules. The section
     stains darker throughout. All 24 dark-only hits in the liver-metastasis group come from
     this section, so its darkness is a staining property. The reviewer excludes the five
     void failures and five named dark fields (MANUAL below). The other eighteen dark fields
     stay in the analysis.

  5. EX-17 is excluded as a whole section. It has 3 analysed fields and the gate drops
     1, which leaves too few for a per-section mean. qc_excluded_samples.csv and an
     _EXCLUDED folder suffix carry that exclusion, so this script skips the section.

The script preserves existing exclusions and never rewrites them. Earlier reviewer decisions
therefore keep their recorded reasons. Re-running is safe because the script skips any field
that is already excluded.

Usage:
    python3 scripts/apply_crop_gate_exclusions.py --dry-run   # show what would change
    python3 scripts/apply_crop_gate_exclusions.py
    bash scripts/rebuild_all_outputs.sh                       # run after applying

Writes: per_sample_out/<subfolder>/<subfolder>_fov_exclude.csv (rows appended per section)
"""

import csv
import glob
import json
import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
DRY = "--dry-run" in sys.argv

MAX_LARGEST_VOID = 0.02      # largest connected near-white component / field area
MIN_TISSUE_FRAC = 0.80       # tissue pixels / field area
MAX_DARK_FRAC = 0.010        # dark pixels / field area

# Fields the reviewer excludes in the one section adjudicated field by field.
MANUAL = {"LM-09": {9, 11, 13, 15, 50, 32, 33, 44, 48, 49}}
SECTION_EXCLUDED = {"EX-17"}

SC = json.load(open("outputs/tables/crop_scores_all.json"))
COH = list(csv.DictReader(open("outputs/tables/cohort_table.csv")))
FLAG = defaultdict(set)
for _r in csv.DictReader(open("outputs/tables/fov_quality_flagged.csv")):
    FLAG[_r["sample_id"]].add(int(float(_r["fov_id"])))


def verdict(s):
    """(should_drop, reason) for one crop's metrics."""
    why = []
    if s["largest_void"] >= MAX_LARGEST_VOID:
        why.append(f"largest void {s['largest_void']:.3f} of field (>= {MAX_LARGEST_VOID})")
    if s["tissue_frac"] <= MIN_TISSUE_FRAC:
        why.append(f"tissue fraction {s['tissue_frac']:.3f} (<= {MIN_TISSUE_FRAC})")
    if s["dark_frac"] >= MAX_DARK_FRAC:
        why.append(f"dark fraction {s['dark_frac']:.4f} (>= {MAX_DARK_FRAC})")
    return bool(why), "; ".join(why)


added_total = skipped_total = files = 0
per_group = defaultdict(int)
for r in COH:
    sid, sub = r["sample_id"], r["subfolder"]
    if sid in SECTION_EXCLUDED:
        continue
    sc = SC.get(sid, {})
    if not sc:
        continue
    path = os.path.join("per_sample_out", sub, f"{sub}_fov_exclude.csv")
    existing = {}
    if os.path.exists(path):
        for x in csv.DictReader(open(path)):
            existing[str(int(float(x["fov_id"])))] = x.get("reason", "")

    mag = {}
    for g in glob.glob(f"per_sample_out/{sub}/fov/*_fov.csv"):
        for x in csv.DictReader(open(g)):
            mag[int(float(x["fov_id"]))] = 1
    candidates = [k for k in mag if k not in FLAG[sid]]

    added = 0
    for k in sorted(candidates):
        key = str(k)
        if key in existing:
            skipped_total += 1
            continue
        if sid in MANUAL:
            if k not in MANUAL[sid]:
                continue
            s = sc.get(key, {})
            reason = ("reviewer-adjudicated 2026-09-02 (crop gate, field-by-field): this "
                      "section stains darker throughout, so its dark fields were reviewed "
                      f"individually; largest void {s.get('largest_void', 0):.3f}, dark "
                      f"fraction {s.get('dark_frac', 0):.4f}")
        else:
            s = sc.get(key)
            if not s:
                continue
            drop, why = verdict(s)
            if not drop:
                continue
            reason = f"reviewer-approved 2026-09-02: crop gate [{why}]"
        existing[key] = reason
        added += 1
        per_group[r["tissue_type"]] += 1

    if added:
        files += 1
        added_total += added
        if not DRY:
            with open(path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["fov_id", "reason"])
                w.writeheader()
                for key in sorted(existing, key=int):
                    w.writerow({"fov_id": key, "reason": existing[key]})
        print(f"  {sid:12s} +{added:3d} excluded  (file now holds {len(existing)})")

print(f"\n  {'DRY RUN -- ' if DRY else ''}{added_total} field(s) newly excluded across "
      f"{files} section(s); {skipped_total} already excluded and left untouched")
for g, n in sorted(per_group.items()):
    print(f"    {g:16s} {n:4d}")
if not DRY:
    print("\n  NEXT: bash scripts/rebuild_all_outputs.sh")
