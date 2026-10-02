# hpc/ — cluster job wrappers

Each wrapper is a Slurm array script. It sources `config/site.env`, runs one Fiji macro from
`macros/` per task, and writes into the project tree. Resources and the Slurm account sit in
each wrapper's `#SBATCH` header.

| Wrapper | Runs | One task per line of | Writes |
|---|---|---|---|
| `run_il15_array.sh` | `il15_ihc_pipeline.ijm` | `imagelist.txt`, `subdirs.txt` | `per_sample_out/<section>/`: slide CSV, field CSV, mask QC, channel figures |
| `run_fov_raw_crops.sh` | `il15_fov_raw_crops.ijm` | `imagelist.txt`, `subdirs.txt` | `per_sample_out/<section>/fov/raw/*_raw.png`, the crops field QC scores |
| `run_fov_quality.sh` | `scripts/pipeline/fov_quality_score.py` | one task | `outputs/tables/fov_quality_*.csv` |
| `run_il15_threshold_scan.sh` | `il15_threshold_scan.ijm` | `inputs_lists/thrscan_samples.txt` | `thrscan/<section>/` |
| `run_dab_threshold_check.sh` | `il15_threshold_scan.ijm` | `inputs_lists/dab_check_inputs.txt`, or `$DAB_LIST` | `outputs/threshold_selection/<DAB_OUT>/` |

`imagelist.txt` holds one source image path per line and `subdirs.txt` the matching output
folder name. They and `inputs_lists/` are not in this repository, because they carry specimen
identifiers and cluster paths. Create your own for your images.

## Running

```bash
sbatch --array=1-74%14 hpc/run_il15_array.sh      # from the project root
```

- Line *N* of `subdirs.txt` is array index *N*.
- Smoke-test one task after any macro edit. A macro error under `xvfb` opens an invisible
  dialog, and the task hangs at 0% CPU until the wall clock.
- Never pass `--headless`. Colour Deconvolution2 builds an AWT dialog and writes nothing
  without a display.
- `IL15_ROOT` overrides the project root, which each wrapper otherwise sets to the cluster
  path. Slurm runs a copy of the wrapper from its spool directory, so the wrapper cannot find
  the project from its own location.
