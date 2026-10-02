#!/usr/bin/env python3
"""Aggregate per-sample IHC results into one cohort table.

Each sample row carries:
  - whole-slide magenta / brown area fraction (single value, from _absorbance.csv)
  - FOV-mean magenta / brown area fraction, the mean over surviving FOVs
    (fov.csv minus any fov_id listed in <sample>_fov_exclude.csv)
Rows are grouped by control_group (brain / liver / liver_met / met) from manifest.csv.

Writes: outputs/tables/cohort_table.csv
"""
import csv, glob, os, re, statistics as st, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from site_paths import field_exclusion_rows
OUT  = os.path.join(ROOT, "outputs","tables","cohort_table.csv")

def sid(s):
    m = re.search(r"(S\d+-\d+)", s); return m.group(1) if m else s

# Map sample_id to (control_group, tissue_type). Key on both the literal manifest id and
# its sid() form so folder names like 'LV-04_3730' (manifest 'LV-04 (3730)')
# resolve. IDs stay distinct under sid(), so the two keys never collide.
meta = {}
with open(os.path.join(ROOT, "manifest.csv"), newline="") as f:
    for r in csv.DictReader(f):
        key = r["sample_id"].strip(); val = (r["control_group"].strip(), r["tissue_type"].strip())
        meta[key] = val
        meta.setdefault(sid(key), val)

# Sample-level QC decisions live in data/decisions/qc_excluded_samples.csv
# (columns sample_id, action=exclude|caveat, reason).
qc = {}
qcf = os.path.join(ROOT, "data", "decisions", "qc_excluded_samples.csv")
if os.path.exists(qcf):
    for r in csv.DictReader(open(qcf)):
        qc[sid(r["sample_id"].strip())] = r.get("action", "").strip()

rows = []
for sdir in sorted(glob.glob(os.path.join(ROOT, "per_sample_out", "*"))):
    if not os.path.isdir(sdir): continue
    sub = os.path.basename(sdir)
    absc = glob.glob(os.path.join(sdir, "*_absorbance.csv"))
    if not absc: continue
    ws = list(csv.DictReader(open(absc[0])))[0]
    # Take the sample ID from the subfolder name, which is built from the manifest.
    # The image column is unreliable because some internal filenames disagree with the
    # manifest ID (e.g. folder EX-03 holds internal file EX-03; folder
    # "LV-04 (3730)").
    sample = sid(sub)
    action = qc.get(sample, "")
    if action == "exclude":
        continue                         # excluded samples leave the table and every stat/figure
    grp, tt = meta.get(sample, ("??", "??"))

    # The FOV mean covers surviving FOVs only.
    fovcsv = glob.glob(os.path.join(sdir, "fov", "*_fov.csv"))
    excl = set()
    for r in field_exclusion_rows(sub):
        v = (r.get("fov_id") or "").strip()
        if v.isdigit(): excl.add(int(v))
    mags, brns, n_tot = [], [], 0
    if fovcsv:
        for r in csv.DictReader(open(fovcsv[0])):
            n_tot += 1
            if int(r["fov_id"]) in excl: continue
            mags.append(float(r["magenta_area_fraction"]))
            brns.append(float(r["brown_area_fraction"]))
    fmean = lambda a: (sum(a)/len(a)) if a else ""
    # Metric 2 is magenta intensity, the mean OD over positive pixels. It is crosstalk-safe.
    mag_pos = None
    try: mag_pos = float(ws.get("mag_pos_mean",""))
    except: mag_pos = None
    # Metric 3 is hematoxylin-normalized intensity, mag_pos_mean / mean hematoxylin OD.
    hema = None
    try: hema = float(ws.get("hema_mean_unfilled",""))
    except: hema = None
    mag_pos_hemanorm = (mag_pos/hema) if (mag_pos is not None and hema and hema>0) else ""
    rows.append({
        "sample_id": sample, "control_group": grp, "tissue_type": tt, "subfolder": sub,
        "ws_magenta_af": float(ws["magenta_area_fraction"]),
        "ws_brown_af":   float(ws["brown_area_fraction"]),
        "ws_mag_pos_mean": mag_pos if mag_pos is not None else "",
        "ws_hema_mean": hema if hema is not None else "",
        "ws_mag_pos_hemanorm": mag_pos_hemanorm,
        "fov_mean_magenta_af": fmean(mags),
        "fov_mean_brown_af":   fmean(brns),
        "fov_sd_magenta_af": (st.pstdev(mags) if len(mags) > 1 else 0.0) if mags else "",
        "n_fov_total": n_tot, "n_fov_used": len(mags), "n_fov_excluded": len(excl),
        "qc_flag": action,
    })

GRP_ORDER = {"brain":0, "liver":1, "liver_met":2, "met":3}
rows.sort(key=lambda r: (GRP_ORDER.get(r["control_group"], 9), r["sample_id"]))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
cols = ["sample_id","control_group","tissue_type","subfolder","ws_magenta_af","ws_brown_af",
        "ws_mag_pos_mean","ws_hema_mean","ws_mag_pos_hemanorm",
        "fov_mean_magenta_af","fov_mean_brown_af","fov_sd_magenta_af",
        "n_fov_total","n_fov_used","n_fov_excluded","qc_flag"]
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

# Print a per-group summary to the console.
print(f"wrote {OUT}  ({len(rows)} samples)")
for g in ["brain","liver","liver_met","met"]:
    gr = [r for r in rows if r["control_group"] == g]
    if not gr: continue
    ws = sorted(r["ws_magenta_af"] for r in gr)
    print(f"  {g:10s} n={len(gr):2d}  whole-slide magenta AF: "
          f"min={min(ws):.5f} med={st.median(ws):.5f} max={max(ws):.5f}")
