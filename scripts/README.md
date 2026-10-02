# scripts/ — workstation analysis

Python that turns the per-section Fiji outputs into the cohort table, statistics and the final
figure. Each script finds the project root from its own location. `site_paths.py` supplies
site paths from `config/site.env` and locates the decision files in `data/decisions/`.

`bash scripts/rebuild_all_outputs.sh` runs the steps marked **R**, in dependency order. Run it
after any change to the per-section results, to `data/decisions/` or to a threshold.

## pipeline/

| Script | Reads | Produces | |
|---|---|---|---|
| `aggregate_cohort.py` | `per_sample_out/*/` slide and field CSVs, `manifest.csv`, `data/decisions/` | `outputs/tables/cohort_table.csv` | R |
| `plot_cohort_stats.py` | `cohort_table.csv` | `outputs/tables/cohort_stats.txt` and two working figures | R |
| `fov_quality_score.py` | native field crops | `outputs/tables/fov_quality_{scores,flagged,by_sample}.csv` (field QC round 1) | |
| `score_all_crops.py` | native field crops | `outputs/tables/crop_scores_all.json` (field QC round 2) | |

`manifest.csv` maps each section to its group. It is not in this repository because it carries
specimen identifiers.

## figures/

| Script | Reads | Produces | |
|---|---|---|---|
| `plot_cohort_bars_pct.py` | `cohort_table.csv` | the final figure, Fig. 3 F and G: `outputs/figures/cohort/fig_cohort_bars_pct_{legend,nolegend}.*`; `outputs/tables/cohort_stats_bars_pct.txt` | R |

## threshold_evidence/

| Script | Reads | Produces | |
|---|---|---|---|
| `summarize_threshold_scan.py` | the Vector Red scan from `hpc/run_il15_threshold_scan.sh` | `threshold_scan_{evidence,summary}.csv` | R |
| `summarize_dab_check.py` | the DAB check from `hpc/run_dab_threshold_check.sh` | `dab_check_summary.{csv,txt}` | |
| `dab_threshold_sensitivity.py` | the DAB full-cohort sweep, `cohort_table.csv` | `dab_sensitivity.{csv,txt}` | |

## qc_record/

The field-QC exclusion rules, as applied. See `qc_record/README.md`.

## Top level

| File | Does |
|---|---|
| `rebuild_all_outputs.sh` | runs the steps marked R |
| `site_paths.py` | project root, site paths, and the field-exclusion helpers `field_exclusions()` and `field_exclusion_rows()` |
