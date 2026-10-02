# Materials and Methods — IL-15 / IL-15Rα duet IHC image quantification

Publication-ready text for the Materials & Methods section, plus the parameter tables a
reviewer or a reproducing lab would need. **Last updated 2026-09-10**, describing the
current analysis (measured stain vectors, chromaticity tissue gate, FOV-mean primary
endpoint, reduced field size for thin liver-metastasis cores, and the two-round field
curation with its 2026-09-10 correction).

> **Fill-in markers.** Items marked `⟦TO CONFIRM⟧` are wet-lab or software-version details
> that were not derivable from the analysis code and must be completed from lab records
> before submission. Nothing in this document is inferred or assumed — if a value was not
> verifiable from the pipeline or its outputs, it is marked rather than estimated.

---

## 1. Specimens and cohort

Formalin-fixed paraffin-embedded (FFPE) **human** tissue sections were analysed
(n = 74 sections collected; n = **57** after quality control, see §8). The cohort comprised
four groups:

| Group | Label used in figures | n collected | n analysed |
|---|---|---|---|
| Normal brain | Non-neoplastic brain | 10 | 10 |
| Breast-cancer brain metastasis | BCBM | 38 | 26 |
| Normal liver | Non-neoplastic liver | 10 | 10 |
| Liver metastasis | Liver metastasis | 16 | 11 |

No normal brain or normal liver section was excluded; all 17 exclusions fall in the two
tumour groups, where specimen quality (fragmentation, needle-core size, slide defects) is
the limiting factor.

Normal brain served as the **negative control** for the analysis: it defines the apparent
IL-15Rα signal produced by colour-deconvolution crosstalk in tissue lacking the target
(§9), rather than functioning as a biological comparator.

**Species corrected 2026-09-02.** This section previously read "mouse tissue sections". The
cohort is **human**: the staining protocol is headed "Human Normal and Metastasis Samples" and
the specimen identifiers (BR-01, BM-12, …) are human surgical-pathology accessions. The
project directory name (`7_mouse_IL-15Ra_Liver_Brain`) is a legacy of an earlier mouse
experiment in the same folder tree and does not describe this cohort.

**Ethics.** Archival human tissue was obtained under a protocol approved by the University of
California, Irvine Institutional Review Board, **UCI 17-05**. The study was conducted in
accordance with the Declaration of Helsinki.

`⟦TO CONFIRM⟧` consent status — whether written informed consent was obtained, or the
requirement waived for retrospective use of existing de-identified material. One of the two
must be stated; journals will not clear the Methods without it. Also confirm that UCI 17-05 is
the protocol covering **this secondary use**, not only the original collection, and check the
number against the approval letter (UCI protocols often appear in a longer form such as
`HS# 2017-XXXX`); patient demographics as reportable; section
thickness; whether sections are serial or independent; specimens per group.

## 2. Immunohistochemistry

Sequential (duet) chromogenic immunohistochemistry was performed for IL-15 and IL-15Rα on
a single section, with haematoxylin counterstain:

- **IL-15** — visualised with **3,3′-diaminobenzidine (DAB)**, brown.
- **IL-15Rα** — visualised with **Vector Red / alkaline-phosphatase substrate**, magenta.
- **Counterstain** — haematoxylin, blue.
- **Detection** — ImmPRESS polymer detection system.

`⟦TO CONFIRM⟧` primary antibody clones/catalogue numbers, host species, dilutions and
incubation times; antigen-retrieval buffer, pH, temperature and duration; blocking
reagents; ImmPRESS kit catalogue numbers; DAB and Vector Red development times; and
whether single-stain controls were run. (Single-stain slides for haematoxylin, DAB and
Vector Red alone were run and imaged; the stain vectors in §4 were measured from them on
2026-08-07, replacing the earlier estimated vectors — see §11.)

## 3. Whole-slide imaging

Sections were imaged on a **Keyence** slide scanner at **10× magnification** with
automated stitching, and exported as tiled BigTIFF. Image resolution was
**0.75488 µm per pixel** (confirmed 2026-08-06 against the calibration reported by the Keyence acquisition software; the value used for every physical measurement below). Note the exported TIFF headers carry only placeholder resolution tags (72/96 dpi), so the scale must be taken from the instrument, not the file. Sections
ranged from 0.39 to 352.6 mm² of detected tissue (median 107.3 mm²).

`⟦TO CONFIRM⟧` exact Keyence model (e.g. BZ-X series), objective NA, illumination and
exposure settings, and whether shading correction / white balance was applied at
acquisition.

## 4. Colour deconvolution

Images were processed in **Fiji / ImageJ** using the **Colour Deconvolution2** plugin
(G. Landini) which implements the Ruifrok–Johnston method. Analysis was automated as an
ImageJ macro (`macros/il15_ihc_pipeline.ijm`) and executed non-interactively on a Linux
compute cluster under a virtual X display (`xvfb-run`).

Images were opened via the Bio-Formats importer with file grouping disabled, and flattened
to 24-bit RGB (Colour Deconvolution2 requires RGB input).

Stain vectors (normalised RGB optical-density unit vectors) were determined from
representative regions of the study material:

| Stain | R | G | B |
|---|---|---|---|
| Haematoxylin | 0.61190 | 0.65217 | 0.44749 |
| IL-15 (DAB) | 0.28420 | 0.50621 | 0.81424 |
| IL-15Rα (Vector Red) | 0.18300 | 0.85109 | 0.49209 |

Vectors were **measured from single-stain control slides** (haematoxylin alone, DAB alone,
Vector Red alone) imaged under identical conditions, rather than estimated from duet-stained
material.

**Limitation (important).** The DAB and Vector Red vectors are near-collinear
(dot product = **0.884**). Under near-collinearity the deconvolution cannot cleanly
apportion absorbance between the two chromogens, so *per-pixel magnitude* in either
channel is unreliable. The analysis therefore uses **thresholded positive-area fraction**
as its endpoint and does **not** interpret mean absorbance / staining intensity (§9).

Background intensity (I₀) was measured per image from the unstained glass rather than
assumed, and optical density computed as OD = −log₁₀(I / I₀) per channel.

## 5. Tissue segmentation

A tissue mask was built per whole-slide image, on an 4×-downsampled copy
(3.02 µm per mask pixel) for tractability, as follows:

1. **Total-OD map** — the red, green and blue channel OD images of the *original RGB image*
   were summed. This is computed before, and independently of, colour deconvolution, and is a
   stain-agnostic "how much is here" map, so faintly-stained or haematoxylin-only tissue
   is retained (a single-channel or saturation threshold misses it).
2. **Smoothing** — Gaussian blur, σ = 2 mask pixels, to reconnect sparsely-stained tissue
   so thin structures form connected regions instead of shattering below the size filter.
3. **Tissue threshold** — total OD ≥ **0.10**.
4. **Near-white (glass) rejection — chromaticity gate.** A pixel additionally had to be
   **chromatic**: at least **7 of 255 levels** of separation between its strongest and its
   weakest OD channel (OD 0–0.5 mapped to 0–255, so one level = 0.00196 OD). *Rationale:*
   dim or dusty glass absorbs weakly but **equally** in all three channels, and three small
   per-channel ODs can sum past the 0.10 total-OD threshold, so empty glass is labelled
   tissue — which inflates the denominator and **dilutes** the reported positive-area
   fraction. Glass is therefore near-achromatic (~1 level of separation) whereas even faint
   tissue retains a stain hue (~15 levels). The gate separates them on **colour**, not on
   **magnitude**.

   An earlier magnitude-based gate (maximum single-channel OD ≥ 0.12) was implemented,
   validated, and then **rejected and disabled**: because dim glass and faint tissue have
   the same absorbance magnitude, it could not distinguish them and removed up to 57% of
   genuine tissue on faintly-stained sections. The chromaticity gate replaced it on
   2026-08-06 after validation (recovering BM-01 from 26 to 101 mm² and BM-23 from
   227 to 272 mm² of tissue, both confirmed genuine on the QC overlays, with clean liver
   unchanged) and — critically — **the normal-brain negative-control floor was unchanged at
   0.03% of tissue area**, confirming the gate does not manufacture signal.
5. **Speckle removal** — morphological opening (1 iteration).
6. **Size filter** — connected components < **0.005 mm²** discarded. Sections contained
   1–136 tissue fragments (median 14), so a permissive filter was required to retain
   genuinely fragmented metastatic tissue.
7. **Hole filling** — interior holes were filled, giving a tissue silhouette. Large-lumen
   carve-out was evaluated and **not** used: excluding large interior lumens fragmented
   sparse metastatic tumour unacceptably. The reported denominator is therefore the filled
   silhouette, which includes interior empty space — vessel lumina, ventricles and sinusoids
   count toward tissue area. This inflates the denominator and therefore **understates** every
   positive-area fraction slightly; the bias is conservative and applies identically to all
   groups. Note that the `hole_fraction` column in `*_absorbance.csv` is **identically zero
   for every section** and does not quantify this: it measures filled-minus-unfilled, and with
   the large-lumen carve-out disabled the unfilled mask *is* the filled mask. An earlier
   version of this section wrongly cited that column as making the contribution auditable
   (corrected 2026-09-09).
8. **Edge erosion** — the mask was eroded by **100 µm** to exclude edge/knife artefacts
   and section-margin staining. Erosion was **relaxed progressively** (to 0 µm where
   necessary) for thin or small tissue that erosion would otherwise remove entirely;
   the erosion actually applied is recorded per image (range 0–100 µm across this cohort).

All reported quantities were measured **through the eroded mask**.

## 6. Primary endpoint: IL-15Rα-positive area fraction

A pixel within the tissue mask was scored **IL-15Rα-positive** if its Vector Red
(magenta) channel optical density was ≥ **0.10**. The primary endpoint is the
**positive-area fraction**:

> IL-15Rα-positive area fraction = (IL-15Rα-positive pixels) / (all tissue-mask pixels)

expressed as a percentage of tissue area. This quantity is computed two ways: as a
**per-sample mean over sampled fields of view** (the primary endpoint — see §7) and over the
**whole slide** (reported alongside it). The two agree in direction throughout; where they
differ in statistical significance for a given contrast, both are reported.

The threshold was set at 0.10 from a **threshold scan** (0.02–0.40 in 0.02 steps) run over
22 sections — all 10 normal brains plus 12 signal-bearing sections — and evaluated against the
normal-brain negative control, which contains no IL-15Rα and therefore reads only residual
crosstalk. Separation between normal liver and the brain floor is **flat from 0.10 to 0.15**
(3118× vs ~3100×) while the measured signal doubles over the same interval, so the higher
threshold discarded genuine faint signal without improving specificity; below 0.08 the ratio
falls as background enters. At 0.10 the brain floor is **0.0036% of tissue area (36 ppm)**.
The choice was confirmed visually against the sections: the pixels gained between 0.15 and
0.10 fall on stained parenchyma, not on glass or vessel lumen. The threshold was then frozen
and applied identically to every sample.

An earlier value of 0.15 was derived before the stain vectors were measured (§4). The
estimated vectors then in use were more collinear (DAB·Vector Red = 0.947 vs 0.884), so DAB
absorbance was partly mis-apportioned into the Vector Red channel; the negative-control floor
was **0.03%** under those vectors versus **0.0013%** with the measured set, a ~23-fold
reduction in tissue that expresses no IL-15Rα. The 0.15 threshold had been compensating for
that crosstalk. Values reported here are consequently lower than, and not comparable to,
any figure produced before 2026-08-07.

**IL-15 (DAB).** IL-15-positive area was computed the same way, at DAB OD ≥ **0.30**. Brain
cannot anchor this threshold, because it is the highest IL-15 group, so the floor comes from
the stain controls. At 0.30 the four non-DAB controls (secondary-only brain and liver,
haematoxylin-only and Vector Red-only slides) show at most **0.63%** DAB-positive tissue,
against IL-15 group medians of 25–65% on the 21 scan sections, a ≥ 40-fold separation. The
group comparison does not depend on the exact cut: across all 57 sections, every IL-15
pairwise contrast kept its significance and direction at every threshold from 0.24 to 0.40
(whole-slide endpoint; criterion fixed before the scan). The order of the two closest groups,
BCBM and liver metastasis, reverses within that range; their difference is non-significant
at every threshold. Evidence: `data/threshold_evidence/dab_check_summary.csv` and
`dab_sensitivity.csv`.

## 7. Field-of-view (FOV) sampling

In addition to the whole-slide measurement, each section was sampled with a
**systematic uniform-random grid** of square fields of view (stereology-standard design,
which has better variance properties than independent random placement):

| Parameter | Value |
|---|---|
| FOV size | 662 px = **500 µm** square = **0.25 mm²** for 50 of the 57 analysed sections; **250 µm** square = **0.0625 mm²** (331 px) for 7 thin liver-metastasis cores (see below) |
| Grid pitch | **1250 µm** default; **625 µm** for 1 small liver-metastasis punch and **312 µm** for the 7 thin cores (FOVs are spaced sample points, not a contiguous tiling) |
| Grid origin | single random offset per image, from a fixed per-image seed (reproducible) |
| Edge erosion | **100 µm** default; **40 µm** / **20 µm** for the reduced-geometry samples above |
| Inclusion | FOV retained only if **≥ 99.9% inside the eroded tissue mask** |
| Minimum | if < 3 FOVs qualified, additional placements were sought on a half-FOV step, and finally a single largest-inscribed FOV; the placement route is recorded per FOV |

**Reduced FOV size for thin cores.** Seven liver-metastasis specimens are thin needle/punch
cores. Under the 500 µm field and the requirement that a field lie entirely within the eroded
tissue mask, each yielded only 0–1 usable field, because a thin core is largely edge and most
grid positions straddle the tissue boundary. For these seven sections only, the field was
reduced to 250 µm (pitch 312 µm, erosion 20 µm), raising their combined yield from 4 to 149
fields. This is legitimate because the endpoint is an **area fraction** — positive area divided
by field area — which is dimensionless and therefore scale-free: a 250 µm field estimates the
same quantity as a 500 µm field. What changes is per-FOV variance (smaller fields are
individually noisier) and the number of fields sampled, which largely offset in the per-sample
mean; no directional bias is introduced. Per-sample field sizes and pitches are recorded in
`data/decisions/fov_sampling_overrides.csv`.

For each FOV the same positive-area fraction was computed. The **primary endpoint is the
per-sample FOV-mean**: the unweighted mean of that sample's qualifying fields. It is unweighted
because 99.7% of qualifying fields are ≥ 99.9% tissue, making area-weighting a no-op while
adding an assumption. The whole-slide measurement (§6) is reported alongside it as
corroboration. Because section size varies, the number of fields per sample varies widely
(3–187, median 57), so **per-sample FOV-means are not equally precise**; four sections rest
on 3–4 fields (three liver-metastasis cores and one BCBM section, BM-10) and are
identified as such in the results.

## 8. Quality control

Every section's tissue mask was reviewed against a per-sample QC overlay showing the
detected tissue boundary before erosion, the eroded mask actually measured, and the
original image. Sections were assessed for mask under-/over-segmentation, tissue
fragmentation, tears and folds. Outcomes:

- **17 sections excluded** (n = 74 → **57 analysed**: 10 non-neoplastic
  brain, 10 non-neoplastic liver, 11 liver metastasis, 26 BCBM). Grounds, recorded per sample, were: a tissue mask
  outlining empty glass with tissue too sparse to detect reliably (EX-04); heavy
  fragmentation (e.g. EX-05, 124 fragments over 74.0 mm²; EX-01, 1.94 fragments
  per mm² with only 4 usable fields); a majority of fields artefact-flagged (EX-02,
  71%; EX-10, 6 of 7); slide-preparation defects such as bubbles (EX-13); and one
  specimen that is **physically unmeasurable** — the needle core of EX-16 is only
  ~222 µm across at its widest point, so no field of either size fits inside it; and one
  section (EX-03) yielding only a single usable field in 1.05 mm² of tissue, which cannot
  support a per-sample mean.
- All exclusions were **reviewer-approved before being applied**; none were made
  automatically. Excluded sections are retained on disk in folders suffixed `_EXCLUDED`
  so the exclusion is visible without consulting the table.
- No section is retained "with a caveat": sections previously carried that way
  (EX-07) were subsequently excluded outright by reviewer decision.
- Automated tissue-**fold** detection was implemented and abandoned (it could not be made
  to run reliably) and contributed to no reported result; folds were instead assessed
  visually on the QC overlays.
- **Field-level curation, round 1 (2026-08-12).** Every measured field was scored at native
  resolution for tears (largest near-white connected component > 6% of the field with
  solidity < 0.75, or total void > 20%, or any single void > 15%), folds (elongated dark band
  > 4% of the field, axis ratio > 3.0) and dried tissue (> 25% very dark and desaturated).
  **304 fields** were excluded on this basis.
- **Field-level curation, round 2 (2026-09-02).** While assembling the supplementary
  deconvolution montage, fields that were ragged, torn at the section margin, or reduced to
  strands were found among the *retained* set: they sit under all three of the round-1 tear
  criteria. Three examples that passed round 1 — BM-16 field 218 (total void 0.204, largest
  single void 0.057), LV-05 field 75 (0.190 / 0.111) and LM-08 field 42 (0.104 / 0.064).
  A second scorer was therefore run on each field's raw 662 px crop, measuring the **largest
  single connected near-white component** (excluded at ≥ 2% of the field), the very-dark
  fraction (≥ 1%) and the tissue fraction (≤ 80%). The discriminating quantity is the largest
  *single* void, not total void: ragged fields score 0.047–0.111 and solid fields 0.001–0.006,
  an ~8-fold separation, whereas a total-void rule would have discarded steatosis and liver
  sinusoids and so preferentially stripped the liver groups. **456 fields** were excluded on
  this basis, across 52 sections.
- **Both rounds were adjudicated visually, on the image**, not on measured values — acting
  only on flagged fields that read high would bias every affected sample downward. Automated
  scores were used solely to **rank fields for review**. One section (LM-09) stains darker
  throughout, so its dark-flagged fields were adjudicated individually rather than by rule:
  10 of its 28 flagged fields were excluded and 18 retained. A proposed refinement keying the
  dark criterion on the *shape* of the dark region rather than its area — which would have
  retained nuclei-dense tumour fields — was reviewed and **not adopted**.
- **Effect on results: none of the reported contrasts changed.** All six pairwise verdicts,
  and the Kruskal–Wallis result, are identical before and after round 2; group medians moved
  by less than 0.02 percentage points. One section, BM-13, fell from 19 to 8 fields and was
  retained by reviewer decision; it is the only section whose estimate rests on materially
  less data than before.
- **767 field exclusions** are recorded in total (2 from an earlier individual review, 304
  from round 1, 461 under the round-2 criteria), each written into its section's
  `*_fov_exclude.csv` with the metrics that produced it. Eighteen by-eye rejections made before
  the 2026-08-10 pipeline re-run had been stored as field numbers on the superseded grid; image
  matching against the archived pre-re-run montages showed they no longer denoted the fields
  that were judged, so they were removed and each affected field re-assessed under the current
  criteria (13 restored, 5 excluded under round 2). No reported contrast changed. The complete
  record, one row per excluded field with its reason, is `data/decisions/field_exclusions/`.
- Per-FOV channel montages used for field-level curation are **verified to match the FOV
  grid that produced them**, by content comparison, since re-running the
  pipeline with a changed tissue mask or field size renumbers the grid and would otherwise
  leave curation images keyed to superseded field identities.

QC decisions and their reasons are recorded in `data/decisions/qc_excluded_samples.csv`.

## 9. What is *not* claimed

Stated explicitly because it follows from the stain chemistry, not from the data analysis:

1. **Staining intensity is not interpreted.** Mean magenta OD over positive pixels did not
   differ between groups (all pairwise comparisons non-significant), which is the expected
   consequence of near-collinear DAB/Vector Red vectors (§4). Only *how much area* is
   positive is interpretable, not *how strongly* it stains.
2. **IL-15 (DAB) is compared on positive area only**, for the same reason. Its threshold
   rests on the stain controls and a full-cohort threshold sweep (§6), not on brain, which
   is the highest IL-15 group. No ordering of BCBM against liver metastasis is claimed.
3. **Normal brain is a technical floor, not a biological zero.** Its apparent IL-15Rα-positive
   area (median 0.0015% of tissue on the FOV-mean endpoint; 0.0036% whole-slide) is
   deconvolution crosstalk, and sets the detection floor against which
   other groups should be read.
4. **Area fraction is not cell counting.** The endpoint is the fraction of tissue *area*
   that is positive; it does not resolve individual cells or distinguish membranous from
   cytoplasmic localisation.

## 10. Statistics

Analyses used **Python 3.13.12** with **SciPy 1.17.1** and **NumPy 2.4.3**
(figures: Matplotlib 3.10.8).

Because the data span several orders of magnitude and are strongly right-skewed,
**non-parametric** tests were used throughout:

- **Omnibus:** Kruskal–Wallis across the four groups.
- **Pairwise:** two-sided Mann–Whitney *U* for all six group pairs.
- **Multiplicity:** Benjamini–Hochberg false-discovery-rate correction across the six
  pairwise tests; adjusted *p* values are reported.
- Significance markers: `*` p < 0.05, `**` p < 0.01, `***` p < 0.001 (adjusted).
- Figures show the group **median** with the **interquartile range**, and all individual
  section values overlaid. Median/IQR — not mean/SEM — because the reported tests are
  rank-based, so the figure summarises the data the same way the statistics do; on a
  log axis mean ± SEM is also visually asymmetric and can imply values below zero.

Both the FOV-mean (primary) and whole-slide endpoints were tested. They agree on five of
the six pairwise contrasts. They differ on **BCBM vs liver metastasis**, which is
significant on the whole-slide endpoint (BH-adjusted p = 0.032) but **not significant** on
the FOV-mean endpoint (p = 0.100). Both are reported; the FOV-mean result governs, being
the primary endpoint. Group medians on the FOV-mean endpoint are: non-neoplastic brain
0.0015%, BCBM 0.0709%, liver metastasis 0.6104%,
non-neoplastic liver 5.4362% of tissue area — a
3,548-fold range between non-neoplastic liver and the
negative control.

> **Note for reproduction in GraphPad Prism.** Prism's standard post-hoc after
> Kruskal–Wallis is **Dunn's** test, which controls family-wise error and is more
> conservative than Benjamini–Hochberg FDR. Re-running these data with Dunn's may render the
> marginal comparison (normal liver vs liver metastasis, BH-adjusted p = 0.022 on the
> FOV-mean endpoint)
> non-significant. This is a difference in multiplicity philosophy, not a discrepancy in the
> underlying data.

## 11. Planned refinements (not part of the present analysis)

Recorded so that the current numbers are not mistaken for the final protocol.

**Completed since this section was first written:**

1. ~~Single-stain slides to measure stain vectors empirically~~ — **done** (2026-08-07).
   Haematoxylin-, DAB- and Vector Red-only slides were imaged under identical conditions and
   the vectors in §4 measured from them. Candidate vector sets were compared by the
   **condition number of the 3×3 stain matrix** rather than by pairwise dot products, since
   conditioning is what governs whether the unmixing is numerically stable. The
   near-collinearity limitation (§4) is reduced but not eliminated, so intensity is still
   not interpreted.
2. ~~Threshold validation~~ — **done**. Secondary-only (no primary antibody) control slides
   were imaged for brain and liver and used to confirm that the magenta threshold sits
   above the non-specific/crosstalk floor.

**Still outstanding:**

3. ~~FOV-level curation~~ — **done**. Two review rounds (2026-08-12 and 2026-09-02) plus the
   2026-09-10 correction of 18 rejections recorded on a superseded grid; see §8. Every number
   in this document is post-curation.
4. **Background subtraction using the secondary-only controls** is under consideration as a
   sensitivity analysis. It can produce small negative values on true-negative sections;
   that is expected for a subtractive correction and is not by itself grounds to reject it,
   but it has not been applied to any reported number.

## 12. Data and code availability

Code: <https://github.com/Dlawsonlab/IL-15-IL-15Ra_FFPE_IHC_Analysis_Pipeline>. Paths below are
those of that repository.

| Item | Location |
|---|---|
| Whole-slide analysis macro (segmentation, deconvolution, quantification) | `macros/il15_ihc_pipeline.ijm` |
| Cluster job wrapper for the cohort run | `hpc/run_il15_array.sh` |
| Field QC: native-resolution crops, damage scoring, exclusion rules | `macros/il15_fov_raw_crops.ijm`; `scripts/pipeline/fov_quality_score.py`, `scripts/pipeline/score_all_crops.py`; `scripts/qc_record/` |
| Cohort aggregation and statistics | `scripts/pipeline/aggregate_cohort.py`, `scripts/pipeline/plot_cohort_stats.py` |
| Final figure, panels F and G | `scripts/figures/plot_cohort_bars_pct.py` |
| Threshold selection: Vector Red scan, DAB check and full-cohort sweep | `hpc/run_il15_threshold_scan.sh`, `hpc/run_dab_threshold_check.sh`, `scripts/threshold_evidence/` |
| Decisions read as input: section and field exclusions, sampling overrides, stain vectors | `data/decisions/` |
| Reported tables and threshold evidence | `data/results/`, `data/threshold_evidence/` |
| Software environment: Python lock, Fiji versions, Apptainer recipe and image hash | `environment/` |
| Site paths (Fiji, storage) | `config/site.env` |
| Regenerate the tables and figure from the per-section outputs | `bash scripts/rebuild_all_outputs.sh` |

Source images are not included.

**Software versions (recorded 2026-09-09).** The analysis ran under **ImageJ 1.54p**
(`ij-1.54p.jar`) within Fiji / ImageJ2 **2.16.0**, installed at
`/pub/evaz/fiji/Fiji.app` on the analysis cluster. Colour deconvolution used **Colour
Deconvolution2** (G. Landini), whose JAR carries no manifest and no internal version string —
its only stable identifier is its content hash:

| | |
|---|---|
| Plugin | `colour_deconvolution2.jar` |
| SHA-256 | `1690e891f5aef7bbe116806f425539c6b00d38e5d9bd7fd11e962026886e1002` |
| Verified identical | cluster `Fiji.app/plugins/` and the archived copy at `macros/colour_deconvolution2.jar` |

Cite the hash rather than a version number; "Colour Deconvolution2" alone does not pin a build.

**Version difference, and why it does not affect the result.** Stain-vector selection was
performed separately in **ImageJ 1.54t**, a later build than the 1.54p used for the cohort
analysis. This does not propagate into the results: the deconvolution arithmetic is implemented
in the plugin, not in ImageJ core, and the plugin binary is byte-identical between the two
(same SHA-256 above); the vectors themselves enter the analysis as nine literal constants in
`macros/il15_ihc_pipeline.ijm`, not as anything recomputed at run time. ImageJ core here
performs image I/O and array arithmetic. This has been reasoned about and the binaries checked,
but **not** tested by running one section through both builds and comparing outputs; do that if
a reviewer presses on it.

---

## Appendix — complete parameter table

Every numeric choice in the pipeline, for direct reproduction.

| Parameter | Value | Meaning |
|---|---|---|
| `PIXEL_UM` | 0.75488 | µm per pixel at 10× |
| `MAGENTA_THRESHOLD` | 0.10 | OD cut for IL-15Rα positivity (**primary endpoint**) |
| `BROWN_THRESHOLD` | 0.30 | OD cut for DAB positivity (reported, not interpreted) |
| `TISSUE_OD_THRESHOLD` | 0.10 | min summed OD for tissue |
| `TISSUE_MAXOD_THRESHOLD` | **0 (disabled)** | superseded 2026-08-06 by the chromaticity gate (§5); retained in the macro set to 0 |
| `TISSUE_CHROMA_LEVELS` | 7 | chromaticity gate: min levels between strongest and weakest OD channel |
| `MASK_DOWNSAMPLE` | 4 | mask built at 1/4 scale (3.02 µm per mask pixel) |
| `MASK_SMOOTH_PX` | 2 | Gaussian σ on the OD-sum map, in mask pixels |
| `MIN_FRAGMENT_MM2` | 0.005 | smallest retained tissue fragment |
| `EDGE_ERODE_UM` | 100 | mask erosion (relaxed adaptively; 0–100 µm applied) |
| `EXCLUDE_LARGE_LUMENS` | 0 (off) | large interior lumens are **not** carved out |
| `HOLE_FILL_MAXMM2` | 0.03 | (inactive while lumen carve-out is off) |
| `FOV_SIZE_PX` | 662 (331 for 7 thin cores) | 500 µm square = 0.25 mm²; 250 µm = 0.0625 mm² where overridden |
| `FOV_GRID_UM` | 1250 (625 / 312 where overridden) | systematic grid pitch; per-sample overrides in `data/decisions/fov_sampling_overrides.csv` |
| `FOV_REQUIRE_FULL` | 1 | FOV must be ≥ 99.9% inside the eroded mask |
| `MIN_FOV_PER_SAMPLE` | 3 | top-up target for small sections |
| `VECTOR_DdotM` | 0.884 | DAB · Vector Red collinearity (recorded per image) |
| `MAX_LARGEST_VOID` | 0.02 | round-2 field gate: largest single near-white component / field |
| `MAX_DARK_FRAC` | 0.010 | round-2 field gate: very-dark pixels / field |
| `MIN_TISSUE_FRAC` | 0.80 | round-2 field gate: minimum tissue fraction |
