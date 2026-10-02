# data/

Three kinds of file, one folder each.

| Folder | What | Regenerable? |
|---|---|---|
| `decisions/` | Inputs set by judgment: section exclusions, field exclusions, sampling overrides, stain vectors. The scripts read these. | No. `decisions/README.md` describes each. |
| `results/` | Cohort table, statistics and the round-1 field-QC flags, as reported. | Yes, by the scripts in the table below. |
| `threshold_evidence/` | Summaries of the scans behind the two positivity thresholds. | Yes, by the scripts in the table below, from the cluster scans. |

The scripts read decisions from `data/decisions/`, and write results to `outputs/tables/`
and `outputs/threshold_selection/`. `results/` and `threshold_evidence/` are copies of those
outputs as they stood when the package was built, so a reader can check a rerun against them.

## results/

| File | Contents | Written by |
|---|---|---|
| `cohort_table.csv` | One row per analysed section: group, positive-area fractions (FOV mean and whole slide), field counts | `scripts/pipeline/aggregate_cohort.py` |
| `cohort_stats.txt` | Kruskal–Wallis and pairwise Mann–Whitney U, BH-corrected, for the reported endpoints | `scripts/pipeline/plot_cohort_stats.py` |
| `cohort_stats_bars_pct.txt` | The tests drawn on the final figure, Fig. 3 F–G | `scripts/figures/plot_cohort_bars_pct.py` |
| `fov_quality_flagged.csv` | Fields flagged by the round-1 damage scorer | `scripts/pipeline/fov_quality_score.py` |

## threshold_evidence/

| File | Contents | Written by |
|---|---|---|
| `threshold_scan_summary.csv` | Vector Red group medians (% of tissue) and the liver-to-brain ratio at each of 20 thresholds | `scripts/threshold_evidence/summarize_threshold_scan.py`, from the per-section output of `hpc/run_il15_threshold_scan.sh` |
| `dab_check_summary.csv` | DAB group medians, the control floor and each group's separation from it, per threshold | `scripts/threshold_evidence/summarize_dab_check.py` |
| `dab_sensitivity.csv` | Each IL-15 pairwise contrast (medians, direction, p, BH-adjusted p) at every DAB threshold from 0.24 to 0.40 | `scripts/threshold_evidence/dab_threshold_sensitivity.py` |

`docs/METHODS.md` §6 explains how these fix the two thresholds.
