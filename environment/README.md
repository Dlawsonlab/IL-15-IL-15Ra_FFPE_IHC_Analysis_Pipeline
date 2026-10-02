# environment/ — the software that produced the numbers

Two environments ran the analysis. Fiji ran the macros on the cluster, and Python ran
`scripts/` on a workstation. Both are recorded here, each as a record of what ran and as a
recipe to rebuild it.

## Fiji (cluster)

| File | What |
|---|---|
| `fiji_versions.txt` | ImageJ 1.54p, ImageJ2 2.16.0, Fiji 2.16.1-SNAPSHOT, Bio-Formats 8.1.1, Zulu OpenJDK 1.8.0_452; Colour Deconvolution2 hash, checked against `macros/colour_deconvolution2.jar` |
| `fiji_manifest.tsv` | Path, size and SHA-256 of every file in the Fiji install that ran the cohort |
| `fiji.def` | Apptainer recipe: Rocky Linux 8, Xvfb, and a copy of that install |
| `fiji_image.sha256` | Hash of the image built from `fiji.def` |
| `build_fiji_image.sh` | Slurm job: writes the manifest and versions, builds the image, hashes it |
| `test_fiji_image.sh` | Slurm job: reruns `il15_ihc_pipeline.ijm` on sections inside the image and compares the CSVs with `per_sample_out/` |

The image copies the installed Fiji rather than downloading one, because a fresh download
pulls current update-site jars instead of the ones that ran.

**Verified 2026-09-30.** Two sections were rerun inside the image, a BCBM and a non-neoplastic
liver (array lines 21 and 59). The slide CSV and the per-field CSV of each were byte-identical
to the cohort run. Job 57386791.

**Verified 2026-10-01, with sampling overrides.** Two liver-metastasis cores were rerun the same
way: LM-03 (array line 12; 625 µm grid, 40 µm erosion) and LM-05 (line 32; 250 µm
fields, 312 µm grid, 20 µm erosion). The overrides were read from
`data/decisions/fov_sampling_overrides.csv` with the same `awk` line as `hpc/run_il15_array.sh`.
Both CSVs of each section were byte-identical to the cohort run. Job 57556026.

Run a macro in the image the same way `hpc/run_il15_array.sh` runs it on the host:

```bash
apptainer exec --home "$ROOT" --bind /share/crsp "$SIF" \
    xvfb-run -a /opt/Fiji.app/ImageJ-linux64 --mem=58g --console -port0 \
    -macro "$ROOT/macros/il15_ihc_pipeline.ijm" "<image>|<output dir>/"
```

`test_fiji_image.sh` also sets a private Java prefs directory, which concurrent array tasks
need.

### How the macros were launched

The cohort ran as `xvfb-run -a ImageJ-linux64 --mem=<N>g --console -port0 -macro <macro> "<args>"`,
the command in `hpc/run_il15_array.sh`. Three constraints matter:

- **Never `--headless`.** Colour Deconvolution2 builds an AWT dialog. Under `--headless` it
  throws a HeadlessException and writes nothing; `xvfb-run` supplies the display instead.
- **`-macro` and `-batch` differ only on an error.** The pipeline and channel renders ran with
  `-macro`, the figure, crop and threshold-scan jobs with `-batch`. On a macro error, `-macro`
  opens an invisible dialog and the job hangs at 0% CPU, while `-batch` prints the error.
- **Bio-Formats with `group_files=false`.** Keyence `.tif` exports sit beside `.ktl` and
  `10x_Stitch/` companions, and companion grouping crawls them and hangs the open.

Stain vectors were selected separately in ImageJ 1.54t. That does not affect the results: the
deconvolution arithmetic is in the plugin, whose binary is the same in both, and the vectors
enter the pipeline as nine literal constants in `il15_ihc_pipeline.ijm`.

The image itself (about 600 MB) lives outside the project at
`/pub/evaz/containers/il15_fiji_20260930.sif`. It is too large for git, so deposit it in an
archive such as Zenodo and cite the DOI with the hash in `fiji_image.sha256`.

## Python (workstation)

| File | What |
|---|---|
| `python.yaml` | Conda environment, locked to the workstation that ran `scripts/` (macOS arm64): Python 3.13.12, NumPy 2.4.3, SciPy 1.17.1, Matplotlib 3.10.8, Pillow 12.1.1, openpyxl 3.1.5, tifffile 2026.5.2, imagecodecs 2026.3.6, and everything they pull in |

```bash
conda env create -f environment/python.yaml
```
