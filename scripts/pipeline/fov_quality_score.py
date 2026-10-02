#!/usr/bin/env python3
"""
Field quality scoring at native resolution. It scores tears, folds and dried or over-dark
tissue.

The script scores native-resolution crops. On the downsampled overview a 662 px field is
only ~17-20 px across, and a void test decided by a dozen pixels cannot tell a genuine tear
from a vessel lumen or a liver sinusoid, so thumbnail-resolution scoring flags tissue types
at very different rates.

Each field's magenta value is recorded for reference and never decides a flag. Regional
variation in staining intensity is normal biology. A field is bad when its tissue is
damaged, whatever its value.

The script scores three failure modes:

  TEAR    large near-white void inside the field. Shape separates a tear from vessels and
          sinusoids. A tear is one big irregular void (low solidity), while sinusoids are
          many small round lumens (high solidity). Two size-only escape clauses also flag
          a field that is mostly void, whatever the void's shape.
  FOLD    doubled tissue, seen as an elongated band markedly darker than the surrounding
          tissue. Elongation of the dark component separates it from dense nuclear
          staining, which is compact.
  DRIED   mounting failure, where a large fraction of the field is very dark AND
          desaturated. Dried tissue goes muddy grey-brown. Genuine haematoxylin is blue and
          genuine DAB is brown, and both keep chroma.

Thresholds are conservative because the script builds a short, reliable review queue.
Nothing here excludes a field. Exclusion stays a reviewer decision recorded in
<sample>_fov_exclude.csv.

Usage : python3 scripts/pipeline/fov_quality_score.py [sample_id ...]   (no ids = every sample)
Input : per_sample_out/<sub>/fov/raw/<base>_fovNN_raw.png   (hpc/run_fov_raw_crops.sh)
Output: outputs/tables/fov_quality_scores.csv               (every field, every metric)
        outputs/tables/fov_quality_flagged.csv              (fields failing a quality test)
        outputs/tables/fov_quality_by_sample.csv            (per-sample flag rates)
"""
import csv, glob, os, re, sys

import numpy as np
from PIL import Image
from scipy import ndimage

Image.MAX_IMAGE_PIXELS = None
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from site_paths import field_exclusion_rows
os.chdir(ROOT)

# ---- native-resolution thresholds -------------------------------------------------------
WHITE_LVL   = 218    # min(R,G,B) above this = near-white (void: tear, crack, glass, lumen)
DARK_LVL    = 70     # max(R,G,B) below this = very dark
SAT_LOW     = 28     # max(R,G,B)-min(R,G,B) below this = desaturated (dried, not chromogen)

TEAR_AREA   = 0.06   # largest single void must exceed 6% of the field ...
TEAR_SOLID  = 0.75   # ... and be irregular (solidity below this) to count as a tear
VOID_TOTAL  = 0.20   # escape clause; this much total void is damage regardless of shape
VOID_LARGEST= 0.15   # escape clause; a single void this large is damage regardless of shape

# The tear test is a hard AND of two thresholds, so a field just past either one escapes
# it. A field that is mostly hole is damaged whatever shape the hole is, so the two escape
# clauses flag on size alone. A reviewer comparison on BM-01 set their values.

FOLD_DARK   = 0.04   # elongated dark band covering > 4% of the field ...
FOLD_ELONG  = 3.0    # ... with major/minor axis ratio above this

DRIED_FRAC  = 0.25   # > 25% of the field very dark AND desaturated


def solidity(mask):
    """Return the component's area divided by its convex-hull area."""
    ys, xs = np.nonzero(mask)
    if len(xs) < 3:
        return 1.0
    try:
        from scipy.spatial import ConvexHull
        pts = np.column_stack([xs, ys]).astype(float)
        if len(np.unique(pts[:, 0])) < 2 or len(np.unique(pts[:, 1])) < 2:
            return 1.0
        hull = ConvexHull(pts)
        return float(mask.sum()) / max(hull.volume, 1.0)   # in 2-D, .volume is the area
    except Exception:
        # fall back to bounding-box fill fraction
        bb = (xs.max() - xs.min() + 1) * (ys.max() - ys.min() + 1)
        return float(mask.sum()) / max(bb, 1)


def elongation(mask):
    """Return the major/minor axis ratio from the component's second moments."""
    ys, xs = np.nonzero(mask)
    if len(xs) < 3:
        return 1.0
    x = xs - xs.mean(); y = ys - ys.mean()
    cov = np.cov(np.vstack([x, y]))
    if not np.all(np.isfinite(cov)):
        return 1.0
    ev = np.linalg.eigvalsh(cov)
    ev = np.clip(ev, 1e-6, None)
    return float(np.sqrt(ev[-1] / ev[0]))


def score(path):
    a = np.asarray(Image.open(path).convert("RGB")).astype(np.int16)
    mx = a.max(2); mn = a.min(2)
    n = mx.size

    void = mn > WHITE_LVL
    dark = mx < DARK_LVL
    desat = (mx - mn) < SAT_LOW

    void_total = float(void.mean())
    lab, k = ndimage.label(void)
    tear_area, tear_solid, tear_edge, n_void = 0.0, 1.0, 0, k
    if k:
        sizes = ndimage.sum(void, lab, range(1, k + 1))
        big = int(np.argmax(sizes)) + 1
        comp = lab == big
        tear_area = float(sizes.max()) / n
        tear_solid = solidity(comp)
        tear_edge = int(comp[0, :].any() or comp[-1, :].any()
                        or comp[:, 0].any() or comp[:, -1].any())

    dl, dk = ndimage.label(dark)
    fold_area, fold_elong = 0.0, 1.0
    if dk:
        dsizes = ndimage.sum(dark, dl, range(1, dk + 1))
        dbig = int(np.argmax(dsizes)) + 1
        dcomp = dl == dbig
        fold_area = float(dsizes.max()) / n
        fold_elong = elongation(dcomp)

    dried_frac = float((dark & desat).mean())

    is_tear = ((tear_area > TEAR_AREA and tear_solid < TEAR_SOLID)
               or void_total > VOID_TOTAL
               or tear_area > VOID_LARGEST)
    is_fold = fold_area > FOLD_DARK and fold_elong > FOLD_ELONG
    is_dried = dried_frac > DRIED_FRAC

    return dict(void_total=round(void_total, 4), n_voids=n_void,
                tear_area=round(tear_area, 4), tear_solidity=round(tear_solid, 3),
                tear_touches_edge=tear_edge,
                fold_area=round(fold_area, 4), fold_elongation=round(fold_elong, 2),
                dried_frac=round(dried_frac, 4),
                is_tear=int(is_tear), is_fold=int(is_fold), is_dried=int(is_dried))


coh = {r["sample_id"]: r for r in csv.DictReader(open("outputs/tables/cohort_table.csv"))}
# Read each field's magenta value from the sample's own fov.csv, the source of truth.
prev = {}
for _sid, _c in coh.items():
    for _f in glob.glob(f"per_sample_out/{_c['subfolder']}/fov/*_fov.csv"):
        for r in csv.DictReader(open(_f)):
            prev[(_sid, str(int(r["fov_id"])))] = {
                "magenta_pct": round(100 * float(r["magenta_area_fraction"]), 5)}

rows, flagged = [], []
only = set(sys.argv[1:])
for sid, c in sorted(coh.items()):
    if only and sid not in only:
        continue
    sub = c["subfolder"]
    crops = sorted(glob.glob(f"per_sample_out/{sub}/fov/raw/*_fov*_raw.png"))
    if not crops:
        print(f"  {sid}: no native crops yet, skipped")
        continue
    excl = set()
    for r in field_exclusion_rows(sub):
        excl.add(str(int(r["fov_id"])))
    nf = 0
    for p in crops:
        mm = re.search(r"_fov(\d+)_raw", p)
        if not mm:
            continue
        fid = str(int(mm.group(1)))
        s = score(p)
        p0 = prev.get((sid, fid), {})
        rec = dict(sample_id=sid, group=c["control_group"], fov_id=fid,
                   magenta_pct=p0.get("magenta_pct", ""),
                   already_excluded="yes" if fid in excl else "no", **s)
        rec["verdict"] = "+".join([k for k, v in
                                   (("tear", s["is_tear"]), ("fold", s["is_fold"]),
                                    ("dried", s["is_dried"])) if v]) or "ok"
        rows.append(rec)
        if rec["verdict"] != "ok" and fid not in excl:
            flagged.append(rec); nf += 1
    print(f"  {sid:13s} {len(crops):4d} fields scored, {nf:3d} failing a quality test")

if not rows:
    print("\nNo native-resolution crops found. Run hpc/run_fov_raw_crops.sh first.")
    sys.exit(1)

cols = list(rows[0].keys())
os.makedirs("outputs/tables", exist_ok=True)
with open("outputs/tables/fov_quality_scores.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
order = {"dried": 0, "fold": 1, "tear": 2}
flagged.sort(key=lambda r: (order.get(r["verdict"].split("+")[0], 9), -r["tear_area"]))
with open("outputs/tables/fov_quality_flagged.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(flagged)

# ---- per-sample summary: flagged vs total ---------------------------------------------
# A raw flag count conflates "a large sample with a few bad fields" with "a small sample that
# is mostly damaged". The rate says whether to curate fields or reconsider the slide, and
# n_remaining says whether the sample still supports a mean after the flagged ones go.
from collections import Counter, defaultdict
tot_by  = Counter(r["sample_id"] for r in rows)
flag_by = Counter(r["sample_id"] for r in rows if r["verdict"] != "ok"
                  and r["already_excluded"] == "no")
excl_by = Counter(r["sample_id"] for r in rows if r["already_excluded"] == "yes")
mode_by = defaultdict(Counter)
for r in rows:
    if r["verdict"] != "ok" and r["already_excluded"] == "no":
        for m in r["verdict"].split("+"):
            mode_by[r["sample_id"]][m] += 1
grp_by = {r["sample_id"]: r["group"] for r in rows}
summ = []
for sid in sorted(tot_by):
    n, fl = tot_by[sid], flag_by[sid]
    summ.append(dict(sample_id=sid, group=grp_by[sid], n_fov_total=n, n_flagged=fl,
                     pct_flagged=round(100.0 * fl / n, 1) if n else 0.0,
                     n_already_excluded=excl_by[sid],
                     n_remaining_if_flagged_dropped=n - excl_by[sid] - fl,
                     tear=mode_by[sid]["tear"], fold=mode_by[sid]["fold"],
                     dried=mode_by[sid]["dried"]))
summ.sort(key=lambda r: (-r["pct_flagged"], -r["n_flagged"]))
with open("outputs/tables/fov_quality_by_sample.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summ[0].keys())); w.writeheader(); w.writerows(summ)
print("\n  per-sample rates -> outputs/tables/fov_quality_by_sample.csv")
print(f"  {'sample':13s} {'group':9s} {'flagged/total':>13} {'rate':>7} {'left':>6}")
for r in summ[:20]:
    if r["n_flagged"] == 0: break
    print(f"    {r['sample_id']:13s} {r['group']:9s} "
          f"{str(r['n_flagged'])+'/'+str(r['n_fov_total']):>13} {r['pct_flagged']:6.1f}% "
          f"{r['n_remaining_if_flagged_dropped']:6d}")

print(f"\n  {len(rows)} fields scored at native resolution")
print(f"  failing a quality test (and not already excluded): {len(flagged)}")
for kind in ("tear", "fold", "dried"):
    print(f"    {kind:6s}: {sum(1 for r in flagged if kind in r['verdict'])}")
print("  -> outputs/tables/fov_quality_flagged.csv")
