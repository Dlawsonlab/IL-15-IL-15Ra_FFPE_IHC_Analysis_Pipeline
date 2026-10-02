#!/usr/bin/env python3
"""
Score every analysed field crop in the cohort for damage, for reviewer adjudication.

The script computes the crop-gate metrics for every field crop of every section in the
cohort table, so the gate decision rests on the whole cohort.

Metrics per 662 px crop (scored at 331 px):
    tissue_frac    pixels with min(R,G,B) <= 215
    void_frac      pixels with min(R,G,B) >  215  (near-white glass / lumen / tear)
    largest_void   largest connected near-white component, as a fraction of the field
    n_voids        near-white components larger than 0.4% of the field
    dark_frac      pixels with max(R,G,B) < 70    (residue, heavy fold)
    largest_dark   largest CONNECTED dark component, as a fraction of the field
    n_dark_blobs   dark components larger than 0.02% of the field
    dark_elong     major/minor axis ratio of the largest dark component
    cohesion       fraction of tissue in its single largest connected piece

The three dark-connectivity metrics can separate a fold from a nucleus, which a dark
fraction cannot. A fold is one large elongated region, and dense chromatin is thousands of
tiny ones. The reviewer kept the dark rule on the fraction, so these metrics are recorded
for reference and gate nothing.

The run is resumable. Results merge into outputs/tables/crop_scores_all.json, so an
interrupted run picks up where it stopped.

Usage:  python3 scripts/pipeline/score_all_crops.py
Writes: outputs/tables/crop_scores_all.json
"""
import csv, glob, json, os, re, sys, time
import numpy as np
from PIL import Image
from scipy import ndimage

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
os.chdir(ROOT)
CACHE = "outputs/tables/crop_scores_all.json"
WHITE, DARK, MINFRAG = 215, 70, 0.004


def score(path, px=331, _tries=4):
    """Score one crop, retrying transient mount failures.

    The CRSP SMB mount intermittently returns a short or empty read. A crop that decodes fine
    on either side of the attempt can raise UnidentifiedImageError in the middle of a long
    run. A retry with a short backoff clears it, so one bad read does not abort the pass.
    """
    for attempt in range(_tries):
        try:
            im = Image.open(path)
            im.load()
            im = im.convert("RGB")
            break
        except Exception:
            if attempt == _tries - 1:
                raise
            time.sleep(0.4 * (attempt + 1))
    else:  # pragma: no cover
        raise RuntimeError(path)
    if im.width != px:
        im = im.resize((px, px), Image.LANCZOS)
    a = np.asarray(im).astype(np.int16)
    mn, mx = a.min(axis=2), a.max(axis=2)
    white = mn > WHITE
    tissue = ~white
    n = white.size
    wl, nw = ndimage.label(white)
    wsz = np.bincount(wl.ravel())[1:] if nw else np.array([0])
    tl, nt = ndimage.label(tissue)
    tsz = np.bincount(tl.ravel())[1:] if nt else np.array([0])
    dk = mx < DARK
    dl, nd = ndimage.label(dk)
    dsz = np.bincount(dl.ravel())[1:] if nd else np.array([0])
    elong = 1.0
    if dsz.size:
        ys, xs = np.nonzero(dl == int(dsz.argmax()) + 1)
        if len(xs) > 4:
            ev = np.linalg.eigvalsh(np.cov(np.vstack([xs, ys])))
            if ev[0] > 1e-9:
                elong = float((ev[1] / ev[0]) ** 0.5)
    return dict(tissue_frac=round(float(tissue.mean()), 5),
                void_frac=round(float(white.mean()), 5),
                largest_void=round(float(wsz.max() / n) if wsz.size else 0.0, 5),
                n_voids=int((wsz >= MINFRAG * n).sum()),
                dark_frac=round(float(dk.mean()), 5),
                largest_dark=round(float(dsz.max() / n) if dsz.size else 0.0, 6),
                n_dark_blobs=int((dsz >= 0.0002 * n).sum()),
                dark_elong=round(elong, 2),
                cohesion=round(float(tsz.max() / tsz.sum()) if tsz.sum() else 0.0, 5))


cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
# Seed from the montage run's scores so those crops are not rescored.
seed = "outputs/tables/montage_crop_scores.json"
if os.path.exists(seed):
    for k, v in json.load(open(seed)).items():
        cache.setdefault(k, {}).update({kk: vv for kk, vv in v.items()
                                        if kk not in cache.get(k, {})})

coh = list(csv.DictReader(open("outputs/tables/cohort_table.csv")))
t0, done, total = time.time(), 0, 0
for i, r in enumerate(coh, 1):
    sid, sub = r["sample_id"], r["subfolder"]
    sc = cache.setdefault(sid, {})
    crops = glob.glob(f"per_sample_out/{sub}/fov/raw/*_raw.png")
    total += len(crops)
    new = 0
    for p in crops:
        k = str(int(re.search(r"_fov(\d+)_raw", p).group(1)))
        if k in sc and "dark_elong" in sc[k]:
            continue
        sc[k] = score(p)
        new += 1
    done += new
    if new:
        json.dump(cache, open(CACHE, "w"), indent=0)
    print(f"  [{i:2d}/{len(coh)}] {sid:12s} {len(crops):4d} crops, {new:4d} new "
          f"({time.time()-t0:5.0f}s)", flush=True)

json.dump(cache, open(CACHE, "w"), indent=0)
print(f"\n  {sum(len(v) for v in cache.values())} crops scored across {len(cache)} sections "
      f"({done} new this run, {time.time()-t0:.0f}s) -> {CACHE}")
