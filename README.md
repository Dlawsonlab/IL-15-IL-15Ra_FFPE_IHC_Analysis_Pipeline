# IL-15 / IL-15Rα duet IHC quantification

Quantification of chromogenic duet IHC for IL-15 (DAB) and IL-15Rα (Vector Red) in human
non-neoplastic brain and liver, breast-cancer brain metastases and liver metastases.
Fiji/ImageJ on a Slurm cluster, then Python.

**Status.** The analysis is complete: 57 sections analysed, 17 excluded at section QC, 3,399
fields analysed after field QC. The final figure is panels F and G of the paper's Figure 3.
Both positivity thresholds are validated, Vector Red against the brain negative-control floor
and DAB against the stain controls and a full-cohort threshold sweep (`docs/METHODS.md` §6).
Specimens are named by
pseudonymous code (`BR-01`, `BM-07`, `LV-05`, `LM-02`, `EX-17`), not by accession number.

## Where things are

| Path | What |
|---|---|
| `macros/` | The Fiji macros: `il15_ihc_pipeline.ijm`, the whole-slide analysis; `il15_threshold_scan.ijm`, for choosing a positivity threshold; `il15_fov_raw_crops.ijm`, field crops for QC. `colour_deconvolution2.jar` is the exact plugin build used. |
| `hpc/` | Slurm wrappers, one job each. `hpc/README.md` maps wrapper to macro, input list and output. |
| `scripts/` | Workstation Python: `pipeline/` (cohort table, statistics, field-QC scoring), `figures/` (the final figure), `threshold_evidence/`, and `qc_record/` (the field-QC exclusion rules, as applied). `scripts/README.md` lists what each reads and produces. |
| `config/site.env` | Every site path: Fiji, and the lab share as the cluster and the workstation see it. Edit this to run elsewhere. |
| `data/` | `decisions/`: inputs set by judgment, which the scripts read (section and field exclusions, sampling overrides, stain vectors). `results/` and `threshold_evidence/`: the reported tables. `data/README.md` says which is which. |
| `environment/` | Software versions and how the macros were launched; the pinned Python environment; the Apptainer recipe and hash for the Fiji image. |
| `docs/METHODS.md` | The method, with parameter tables. |

Per-section outputs (6.4 GB) and figures are not included. Both regenerate from this code
given the source images, which are not included.

## What the pipeline does

1. **Deconvolution.** Colour Deconvolution2 separates haematoxylin, DAB and Vector Red. The
   stain vectors were measured on single-stain control slides and chosen on the conditioning
   of the unmixing matrix (`data/decisions/stain_vectors_FINAL.md`). Measure your own on
   single-stain slides before applying the pipeline to other material.
2. **Tissue mask.** Summed optical density ≥ 0.10 and a chromaticity gate (≥ 7/255 levels
   between the strongest and weakest channel). Dim glass absorbs weakly but equally in all
   channels, so brightness alone cannot separate it from faint tissue.
3. **Field sampling.** A systematic uniform-random grid of 500 µm fields at 1250 µm pitch,
   each kept only if ≥ 99.9% inside the eroded mask. Eight small liver-metastasis specimens
   use reduced geometry; the endpoint is a dimensionless area fraction, so field size does not
   bias it.
4. **Endpoint.** Positive-area fraction per field, averaged per section: Vector Red OD ≥ 0.10
   for IL-15Rα and DAB OD ≥ 0.30 for IL-15. Non-neoplastic brain lacks IL-15Rα and sets the
   Vector Red floor (0.0036% of tissue at 0.10).
5. **Field QC.** Two rounds of damage scoring at native resolution, adjudicated on the image.
   `docs/METHODS.md` §8 gives the criteria and counts.

Staining intensity is not interpreted. The DAB and Vector Red vectors are near-collinear
(dot product 0.884), so per-pixel absorbance cannot be split reliably between the chromogens;
only positive area is reported.

## Running

```bash
conda env create -f environment/python.yaml       # Python, pinned
sbatch --array=1-74%14 hpc/run_il15_array.sh      # 1. per-section analysis on the cluster
bash scripts/rebuild_all_outputs.sh               # 2. cohort table, statistics, final figure
```

Step 1 reads `imagelist.txt` and `subdirs.txt` (one section per line) from the project root,
and step 2 reads `manifest.csv`, which maps each section to its group. None of the three is in
this repository, because they carry specimen identifiers and cluster paths. Field QC runs
between the two steps: `hpc/run_fov_raw_crops.sh`, the two scorers in `scripts/pipeline/`,
then the rules in `scripts/qc_record/`.
Step 2 runs on a workstation that mounts the same share. Fiji can run from the Apptainer
image in `environment/` instead of a local install; `environment/README.md` has the command.

## Licence

MIT for this project's code (`LICENSE`).

`macros/colour_deconvolution2.jar` is Colour Deconvolution2, © 2004–2020 Gabriel Landini,
licensed under the GNU GPL version 3 or later (<https://www.gnu.org/licenses/gpl-3.0.html>). It
is redistributed unmodified, with its source (`Colour_Deconvolution2.java`) inside the JAR, so
the deconvolution can be rerun with the exact build used here. The JAR has no version string;
identify it by SHA-256 `1690e891f5aef7bbe116806f425539c6b00d38e5d9bd7fd11e962026886e1002`.
The macros call it through ImageJ and do not link it, so the MIT licence covers this project's
own code only.
