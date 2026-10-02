# macros/ — Fiji macros

| Macro | Does | Wrapper |
|---|---|---|
| `il15_ihc_pipeline.ijm` | Whole-slide analysis: tissue mask, deconvolution, field grid, per-field and per-slide measurements. Produced every reported number. | `hpc/run_il15_array.sh` |
| `il15_threshold_scan.ijm` | Positive-area fraction at every threshold in a scan range, for choosing a positivity threshold. | `hpc/run_il15_threshold_scan.sh`, `hpc/run_dab_threshold_check.sh` |
| `il15_fov_raw_crops.ijm` | A native-resolution crop of every field, for field QC. | `hpc/run_fov_raw_crops.sh` |
| `colour_deconvolution2.jar` | The exact Colour Deconvolution2 build used. Licence and hash in the top-level `README.md`. | — |

The stain vectors are nine literal constants (`VEC_H_*`, `VEC_D_*`, `VEC_M_*`) at lines 138–140 of
`il15_ihc_pipeline.ijm`. They were measured on this lab's single-stain slides
(`data/decisions/stain_vectors_FINAL.md`). Measure your own before using the macro on other
material.

## il15_threshold_scan.ijm is a fork of the pipeline

Mask, deconvolution and measurement code is identical in both. They are kept as two files
because each is the code that produced its own outputs. The fork differs in 74 lines:

| | `il15_threshold_scan.ijm` |
|---|---|
| `THRESHOLD_SCAN` | 1 (the pipeline has 0) |
| Scan range | 0.02–0.40 in steps of 0.02; DAB to 0.80 |
| `MAGENTA_THRESHOLD` | 0.15, unused in scan mode; the cohort value 0.10 is in the pipeline |
| Channel | `SCAN_CHANNEL` magenta / dab / both (argument 9) |
| Glass reference | a fixed `R,G,B` accepted as argument 10 |
| Overlay PNGs | every second step, file names tagged by channel |

The macros' comments cite `outputs/tables/fov_sampling_overrides.csv`, which is now
`data/decisions/fov_sampling_overrides.csv`. The wrappers read the new path and pass the
values in as arguments.
