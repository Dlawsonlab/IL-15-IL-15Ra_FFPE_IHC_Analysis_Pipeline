# scripts/qc_record/ — field-QC exclusion rules, as applied

These two scripts turned the field-QC scores into exclusions. On this dataset they ran once,
on the dates below, and their output is `data/decisions/field_exclusions/`. Rerunning them
here would duplicate rows. On new data they are the rules to apply.

| Script | Applied | Reads | Writes |
|---|---|---|---|
| `apply_fov_quality_exclusions.py` | 2026-08-12 | `fov_quality_flagged.csv`, from `scripts/pipeline/fov_quality_score.py` | round-1 field exclusions (native-resolution damage scorer) |
| `apply_crop_gate_exclusions.py` | 2026-09-02 | `crop_scores_all.json`, from `scripts/pipeline/score_all_crops.py` | round-2 field exclusions (crop gate) |

The scripts are kept as they ran. Their comments cite sections of an internal analysis log that
is not part of this repository, and their paths describe the layout of the time:

| Path in these scripts | Current location |
|---|---|
| `scripts/<name>.py` | `scripts/<role>/<name>.py` |
| `outputs/tables/qc_excluded_samples.csv` | `data/decisions/qc_excluded_samples.csv` |
| `per_sample_out/<section>/<section>_fov_exclude.csv` | `data/decisions/field_exclusions/<section>_fov_exclude.csv` |
