# Stain vectors in use — FINAL

These are the vectors the pipeline actually uses. They are **measured from single-stain
control slides** (2026-08-07), not estimated from duet-stained material.

| stain | Colour_ | R | G | B |
|---|---|---|---|---|
| Haematoxylin | 1 | 0.61190 | 0.65217 | 0.44749 |
| DAB (IL-15) | 2 | 0.28420 | 0.50621 | 0.81424 |
| Vector Red (IL-15Rα) | 3 | 0.18300 | 0.85109 | 0.49209 |

Background I0 at the time of picking: 224.76, 232.21, 226.45.
DAB·Magenta dot product: **0.884** (recorded per-image in `*_absorbance.csv` as `vector_DdotM`).

## Why this set

Selected on the **condition number of the 3×3 unmixing matrix**, not on pairwise angles —
deconvolution inverts that matrix, so joint conditioning governs noise amplification:

| vector set | condition number | \|det\| |
|---|---|---|
| old (estimated) | 7.2 | 0.165 |
| **this set (measured)** | **5.0** | **0.199** |
| hybrid (old H + new D,M) | 4.9 | 0.268 |

The hybrid scored marginally better but mixes one legacy estimated vector with two measured
ones, which is not defensible in a methods section for a negligible gain.

## Why it matters

Adopting these reduced the brain negative-control floor from **0.03% to 0.0013%** of tissue
area — a ~23× drop in tissue that contains no IL-15Rα. The old estimated vectors
(DAB·Magenta 0.947) were mis-apportioning DAB absorbance into the Vector Red channel.
**Any value produced before 2026-08-07 is inflated and not comparable.**

## Where these live

- Authoritative: `macros/il15_ihc_pipeline.ijm` (`VEC_H_*`, `VEC_D_*`, `VEC_M_*`)
- Any other macro that deconvolves must carry the same constants (two field-rendering macros
  once drifted, and review images were rendered with the old vectors for five days)
- Machine-readable copy: `data/decisions/stain_vectors_FINAL.csv`
