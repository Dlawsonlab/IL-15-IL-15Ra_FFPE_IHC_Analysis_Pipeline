// il15_threshold_scan.ijm
// Threshold-scan copy of il15_ihc_pipeline.ijm. The mask, deconvolution and measurement
// code is the same; THRESHOLD_SCAN = 1, and the scan settings and arguments 9-10 differ.
// For each slide it writes the positive-area fraction at every threshold in the scan range,
// plus overlay PNGs, and skips the field grid and channel figures.
//
// Run by hpc/run_il15_threshold_scan.sh (the Vector Red scan behind MAGENTA_THRESHOLD) and
// hpc/run_dab_threshold_check.sh (the DAB check of BROWN_THRESHOLD).
//
// Arguments are "|"-separated. Fill every earlier position with a real value, because
// split() collapses empty fields and would shift later arguments left.
//   1 image  2 output dir  3 grid pitch um  4 edge erosion um  5 min tissue fraction
//   6 field size um  7 chroma levels  8 tissue OD  9 channel magenta|dab|both
//   10 fixed glass reference "R,G,B"
//
// ============================================================================
// READ BEFORE USING ANY NUMBER THIS PRODUCES
//
// The DAB and Vector Red stain vectors are near-collinear (D.M = 0.884; condition number of
// the unmixing matrix 5.0). The vectors were measured on single-stain control slides, and
// vectors/00_FINAL_VECTORS_IN_USE.md has the record. Consequences:
//
//   - Mean absorbance over tissue is not a valid endpoint. Low-OD pixels dominate it, and
//     crosstalk dominates there. The CSV carries it for inspection only.
//
//   - Crosstalk depends on signal, so it differs between tissue types. Heavy Vector Red
//     (liver) inflates the DAB channel, and diffuse DAB (brain) inflates the Vector Red
//     channel. Uniform processing does not cancel this.
//
//   - The endpoint is positive AREA fraction at a fixed threshold. It keys on discrete,
//     high-contrast signal and discards the diffuse low-OD regime.
//
//   - Non-neoplastic brain lacks IL-15Ra and serves as the negative control. What it reads
//     is the crosstalk floor. Never tune a threshold to push brain toward zero.
// ============================================================================
//
// Per slide:
//   1. Open the RGB image and set the pixel calibration (the embedded scale is ignored).
//   2. Measure the glass white reference I0 from the slide corners.
//   3. Build the tissue mask on a downsampled copy: summed-OD threshold AND chromaticity
//      gate -> opening -> fragment filter -> hole fill -> edge erosion.
//   4. Colour Deconvolution2 with the fixed vectors, 32-bit absorbance.
//   5. Positive-area fraction for Vector Red and DAB, over the whole section and per field
//      on a systematic uniform random grid. Mean absorbance per channel, inspection only.
//   6. Write the per-slide CSV, the per-field CSV, mask-QC and threshold-QC PNGs and the
//      channel figures.
//
// No rolling-ball background subtraction. It is a local high-pass filter in intensity space,
// breaks Beer-Lambert additivity and would erase cellular-scale signal. OD is referenced to
// measured glass instead. No per-slide vector picking; the vectors are fixed constants.
//
// MODES
//   THRESHOLD_SCAN = 0 : cohort run at MAGENTA_THRESHOLD and BROWN_THRESHOLD.
//   THRESHOLD_SCAN = 1 : area fraction across a threshold range, plus overlay PNGs.
//                        Skips the field grid and channel figures.
//
// Usage, one slide (hpc/run_il15_array.sh does this under Slurm):
//   xvfb-run -a ImageJ-linux64 --console -macro <this macro> '/path/img.tif|/path/outdir'
// Never pass --headless. Colour Deconvolution2 builds an AWT dialog and writes nothing
// without a display. A folder in place of the image processes every image in it; with no
// argument the macro prompts for folders.

requires("1.53a");

// ============================================================
// ==================== USER PARAMETERS =======================
// ============================================================

// ---- Pixel calibration (Keyence 10x; verified on single-FOV and stitched) ----
PIXEL_UM = 0.75488;          // micron per pixel

// ---- MODE ----
THRESHOLD_SCAN = 1;          // 1 = scan (set threshold). 0 = cohort run.
SCAN_MIN  = 0.02;            // Vector Red range 0.02-0.40 in steps of 0.02
SCAN_MAX  = 0.40;
SCAN_STEP = 0.02;
SCAN_CHANNEL = "magenta";    // "magenta" | "dab" | "both" (arg 9)
DAB_SCAN_MAX = 0.80;         // DAB range 0.02-0.80, bracketing the 0.30 cut on both sides

// ---- PRIMARY ENDPOINT ----
MAGENTA_THRESHOLD = 0.15;    // unused in scan mode; the cohort value is in il15_ihc_pipeline.ijm
THRESHOLD_IS_PLACEHOLDER = 0;  // Read twice below. Without this assignment the macro dies on
                             // an undefined variable and hangs invisibly under xvfb.

BROWN_THRESHOLD = 0.30;      // DAB (IL-15) positive-area threshold, checked against the stain
                             // controls and a 0.24-0.40 sweep (outputs/threshold_selection/).

// ---- FOV grid (systematic uniform random sampling) ----
// Fixed spacing: small samples legitimately yield fewer FOVs. n is reported
// per image so you can weight or exclude downstream.
FOV_SIZE_PX      = 662;      // 662 * 0.75488 = 500 um square (0.25 mm^2)
FOV_GRID_UM      = 1250;     // grid pitch, microns
FOV_MIN_TISSUE_FRAC = 0.999; // minimum tissue fraction for a FOV to be kept. 0.999 keeps
                             // only fields entirely inside the eroded mask.
FOV_REQUIRE_FULL = 1;        // FOV must lie ENTIRELY inside the eroded mask;
                             // a straddling FOV reintroduces the edge effect
                             // erosion just removed.
FOV_SEED_BASE    = 20260717; // random start seeded per image (base + filename
                             // hash) so runs are reproducible.
MIN_FOV_PER_SAMPLE = 3;      // if the systematic grid places fewer than this, top up
                             // with non-overlapping full-size FOVs fully inside tissue
                             // (flagged placement=topup) until reached or none fit.
                             // Set 0 to disable top-up (grid + 1-FOV fallback only).

// ---- Edge erosion ----
EDGE_ERODE_UM          = 100;
EDGE_ERODE_MAXFRAC     = 0.40;
EDGE_ERODE_FALLBACK_UM = 40;

// ---- Tissue detection ----
TISSUE_OD_THRESHOLD = 0.10;  // summed-OD tissue threshold. Low enough for faintly stained
                             // tumour; clean brain area is flat across thresholds.
TISSUE_CHROMA_LEVELS = 7;    // chromaticity gate. A pixel is tissue only if its strongest and
                             // weakest OD channels differ by >= 7 8-bit levels (OD 0-0.5 ->
                             // 0-255; 1 level = 0.00196 OD). Glass is achromatic (~1 level) and
                             // faint tissue keeps a stain hue (~15). A magnitude gate cannot
                             // separate them, because dim glass and faint tissue absorb equally.
TISSUE_MAXOD_THRESHOLD = 0;  // max single-channel OD gate; 0 disables it. The chromaticity
                             // gate does this job.
MIN_FRAGMENT_MM2    = 0.005; // smallest fragment kept. Keeps needle cores and scattered tumour
                             // pieces; any debris it admits carries no stain and only pads
                             // the denominator.
MASK_SMOOTH_PX      = 2;     // Gaussian sigma (mask px, ~12 um) on the OD map before
                             // thresholding -> reconnects faint/sparse tissue that
                             // would otherwise shatter into sub-threshold fragments.
HOLE_FILL_MAXMM2    = 0.03;  // (only used when EXCLUDE_LARGE_LUMENS=1) fill interior
                             // holes SMALLER than this; keep larger lumens excluded.
EXCLUDE_LARGE_LUMENS = 0;    // 0: the primary mask is the hole-filled tissue silhouette
                             // (unfilled == filled). 1 carves out large lumens, which gives
                             // erratic masks on fragmented tumour (BM-07) and removes most
                             // of a glandular tumour as lumen (EX-05, 92%).
MASK_DOWNSAMPLE     = 4;     // mask grid 0.75488 * 4 = 3.02 um/px. A coarser grid smears thin
                             // cores below threshold (LM-05 keeps 1 of 7 cores at x8).

// Folds are judged by eye on the mask-QC overlays. Automated fold detection hung under xvfb
// and is not part of this macro; archive/abandoned/FOLD_DETECTION_abandoned.md has the record.

// ---- Background (I0) ----
I0_MANUAL = 0;
I0_R = 240;  I0_G = 240;  I0_B = 240;

// ---- Stain vectors, measured on single-stain control slides (see header) ----
// Order: 1 = Hematoxylin, 2 = DAB (IL-15), 3 = Vector Red (IL-15Ra).
// H.D = 0.868, H.M = 0.887, D.M = 0.884; condition number of the unmixing matrix 5.0.
// vectors/00_FINAL_VECTORS_IN_USE.md has the measurement record.
VEC_H_R = 0.61190;  VEC_H_G = 0.65217;  VEC_H_B = 0.44749;
VEC_D_R = 0.28420;  VEC_D_G = 0.50621;  VEC_D_B = 0.81424;
VEC_M_R = 0.18300;  VEC_M_G = 0.85109;  VEC_M_B = 0.49209;
VECTOR_DdotM = 0.884;        // recorded into the CSV for provenance

// ---- Output ----
SAVE_CHANNEL_TIFS = 0;       // full-res 32-bit absorbance TIFs (large; ImageJ-only LUT).
SAVE_FIGURE_PNGS  = 0;       // full-res Flatten path; needs a display canvas. Keep 0.
SAVE_STAIN_FIGS   = 1;       // downsampled colour PNG per channel + labelled 3-panel montage.
STAIN_FIG_WIDTH   = 1200;    // per-panel width (px) in the stain figures.
SAVE_MASK_QC      = 1;
SAVE_THRESH_QC    = 1;       // keep on -- this is how you see whether the
                             // threshold is catching magenta or brown bleed.
SAVE_FOV_OVERLAY  = 1;       // whole-slide thumbnail with every FOV drawn as a
                             // numbered box (number = fov_id in <base>_fov.csv).
FOV_OVERLAY_WIDTH = 2400;    // gives ~40 px field boxes on a large slide

// ---- LUT display ranges (visualization only) ----
LUT_H_MIN = 0.00;  LUT_H_MAX = 1.20;
LUT_D_MIN = 0.00;  LUT_D_MAX = 1.00;
LUT_M_MIN = 0.00;  LUT_M_MAX = 1.00;

// ---- Stain display colours (white at 0 absorbance -> this at max). Truer to
//      the actual chromogens so it is obvious which panel is which stain. ----
COL_H_R = 60;   COL_H_G = 80;   COL_H_B = 170;   // Hematoxylin  : blue-purple
COL_D_R = 140;  COL_D_G = 80;   COL_D_B = 30;    // DAB / IL-15   : warm brown
COL_M_R = 200;  COL_M_G = 40;   COL_M_B = 110;   // Vector Red / IL-15Ra : magenta-red

// ============================================================
// ======================= DRIVER =============================
// ============================================================

setBatchMode(true);
setOption("BlackBackground", true);
run("Colors...", "foreground=white background=black selection=yellow");

arg = getArgument();
inPath = ""; outDir = "";

if (arg != "") {
    parts = split(arg, "|");
    if (parts.length < 2) exit("Argument must be: <inputPathOrDir>|<outputDir>" +
                               "[|<fovGridUm>|<edgeErodeUm>|<fovMinTissueFrac>]");
    inPath = String.trim(parts[0]);
    outDir = String.trim(parts[1]);
    // ---- optional per-sample sampling overrides ----
    // Small biopsy cores need a tighter grid, less erosion and sometimes a smaller field
    // than whole sections (outputs/tables/fov_sampling_overrides.csv).
    // The ImageJ macro language does NOT short-circuit &&. A guard written as
    //   if (parts.length >= 3 && String.trim(parts[2]) != "")
    // still evaluates parts[2] when only 2 args are passed. The index error opens a modal
    // dialog that is invisible under xvfb, and the job hangs at 0% CPU with no output.
    // Keep these conditions nested.
    if (parts.length >= 3) {
        if (String.trim(parts[2]) != "") {
            FOV_GRID_UM = parseFloat(String.trim(parts[2]));
            print("  [override] FOV_GRID_UM (grid pitch) -> " + FOV_GRID_UM + " um");
        }
    }
    if (parts.length >= 4) {
        if (String.trim(parts[3]) != "") {
            EDGE_ERODE_UM = parseFloat(String.trim(parts[3]));
            print("  [override] EDGE_ERODE_UM -> " + EDGE_ERODE_UM + " um");
        }
    }
    if (parts.length >= 5) {
        if (String.trim(parts[4]) != "") {
            FOV_MIN_TISSUE_FRAC = parseFloat(String.trim(parts[4]));
            print("  [override] FOV_MIN_TISSUE_FRAC -> " + FOV_MIN_TISSUE_FRAC);
        }
    }
    // arg 6 = FOV SIZE in microns. Needle cores can be narrower than a 500 um field, so no
    // 500 um square fits entirely inside the eroded mask and the sample yields 1-3 FOVs.
    // A smaller field recovers real sampling. Positive-AREA-FRACTION is dimensionless, so a
    // 250 um field estimates the same quantity as a 500 um field (only per-FOV variance
    // changes) -- the per-sample mean stays unbiased and cross-group comparison remains valid.
    if (parts.length >= 6) {
        if (String.trim(parts[5]) != "") {
            FOV_SIZE_PX = round(parseFloat(String.trim(parts[5])) / PIXEL_UM);
            print("  [override] FOV size -> " + String.trim(parts[5]) + " um = " + FOV_SIZE_PX + " px");
        }
    }
    // args 7 and 8: chromaticity gate level and total-OD tissue threshold.
    // Callers must pass a real value in every earlier position, because ImageJ's split()
    // collapses consecutive delimiters and blank placeholders shift arguments left.
    if (parts.length >= 7) {
        if (String.trim(parts[6]) != "") {
            TISSUE_CHROMA_LEVELS = parseFloat(String.trim(parts[6]));
            print("  [override] chroma levels -> " + TISSUE_CHROMA_LEVELS);
        }
    }
    if (parts.length >= 8) {
        if (String.trim(parts[7]) != "") {
            TISSUE_OD_THRESHOLD = parseFloat(String.trim(parts[7]));
            print("  [override] tissue OD threshold -> " + TISSUE_OD_THRESHOLD);
        }
    }
    if (parts.length >= 9) {
        if (String.trim(parts[8]) != "") {
            SCAN_CHANNEL = String.trim(parts[8]);
            print("  [override] scan channel -> " + SCAN_CHANNEL);
        }
    }
    // Arg 10: fixed glass reference "R,G,B". For fields that are tissue edge to edge, where
    // the corners are not glass (the single-stain control fields).
    if (parts.length >= 10) {
        if (String.trim(parts[9]) != "") {
            i0p = split(String.trim(parts[9]), ",");
            I0_MANUAL = 1;
            I0_R = parseFloat(i0p[0]); I0_G = parseFloat(i0p[1]); I0_B = parseFloat(i0p[2]);
            print("  [override] fixed I0 -> " + I0_R + ", " + I0_G + ", " + I0_B);
        }
    }
} else {
    inPath = getDirectory("Choose INPUT folder (or cancel to pick a single file)");
    if (inPath == "") inPath = File.openDialog("Choose a single input image");
    if (inPath == "") exit("No input chosen.");
    outDir = getDirectory("Choose OUTPUT folder");
    if (outDir == "") exit("No output folder chosen.");
}
if (!endsWith(outDir, File.separator)) outDir = outDir + File.separator;
File.makeDirectory(outDir);
if (!File.exists(outDir)) exit("Cannot create output folder: " + outDir);

print("");
if (THRESHOLD_SCAN) {
    print("=== THRESHOLD SCAN MODE ===");
    print("Magenta area fraction, " + SCAN_MIN + " to " + SCAN_MAX + " step " + SCAN_STEP + ".");
    print("Run on 2-3 LIVER sections. Compare _scan_thrX.XX.png against the source");
    print("and pick the threshold matching what you call a positive cell.");
    print("Then: set MAGENTA_THRESHOLD, THRESHOLD_IS_PLACEHOLDER=0, THRESHOLD_SCAN=0.");
    print("Do NOT choose the threshold by looking at brain output.");
} else if (THRESHOLD_IS_PLACEHOLDER) {
    print("!!! MAGENTA_THRESHOLD is still flagged PLACEHOLDER.");
    print("!!! Run THRESHOLD_SCAN=1 first. Area fractions below are not meaningful.");
}
print("");

if (File.isDirectory(inPath)) {
    if (!endsWith(inPath, File.separator)) inPath = inPath + File.separator;
    list = getFileList(inPath);
    Array.sort(list);
    n = 0;
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) continue;
        if (!isImageFile(list[i])) continue;
        processSlide(inPath + list[i], outDir);
        run("Close All");
        n++;
    }
    print("Folder complete: " + n + " image(s) -> " + outDir);
} else {
    processSlide(inPath, outDir);
    run("Close All");
    print("Complete: " + inPath);
}

setBatchMode(false);
if (arg != "") eval("script", "java.lang.System.exit(0);");


// ============================================================
// ==================== MAIN PER-SLIDE ========================
// ============================================================

function processSlide(path, outDir) {
    t0 = getTime();
    base = File.getNameWithoutExtension(path);
    print("--- " + base + " ---");

    // ---------- 1. Open ----------
    // Keyence exports are TILED, uncompressed BigTIFFs. ImageJ's native
    // TiffDecoder is strip-based and cannot read tiled TIFFs; a plain open()
    // silently delegates to Bio-Formats, whose import dialog then BLOCKS
    // FOREVER under a virtual display (xvfb). So call Bio-Formats directly,
    // non-interactively (full option string = no dialog).
    // group_files=false is required. The raw Keyence .tif sits next to its .ktl and
    // 10x_Stitch/ companions, and Bio-Formats grouping crawls them and hangs at
    // "Populating OME metadata".
    run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Default " +
        "view=Hyperstack stack_order=XYCZT use_virtual_stack=false group_files=false");
    if (nImages == 0) { print("  FAILED to open: " + path); return; }

    src = getTitle();
    // Bio-Formats returns RGB TIFFs as a 3-channel 8-bit stack, not 24-bit RGB.
    // Flatten to true RGB (Colour Deconvolution2 requires a 24-bit RGB input).
    // Stack to RGB yields a "<title> (RGB)" window -- findDeconWindow() already
    // matches that suffix.
    if (bitDepth() != 24) {
        getDimensions(W0, H0, ch0, sl0, fr0);
        if (ch0 >= 3) {
            run("Stack to RGB");
            rgbTitle = getTitle();
            close(src);
            selectWindow(rgbTitle);
            src = rgbTitle;
        } else {
            run("RGB Color");
            src = getTitle();
        }
    }
    // Overwrite embedded calibration -- do not trust file metadata.
    run("Set Scale...", "distance=1 known=" + PIXEL_UM + " unit=micron");
    getDimensions(W, H, ch, sl, fr);
    print("  " + W + " x " + H + " px  (" + d2s(W*PIXEL_UM/1000,2) + " x " +
          d2s(H*PIXEL_UM/1000,2) + " mm)");

    // ---------- 2. White reference ----------
    i0 = measureI0(src, W, H);
    print("  I0 (R,G,B) = " + d2s(i0[0],1) + ", " + d2s(i0[1],1) + ", " + d2s(i0[2],1));

    // ---------- 3. Tissue masks ----------
    maskInfo = buildMasks(src, W, H, i0, base, outDir);
    tissueArea_mm2 = maskInfo[0];
    filledArea_mm2 = maskInfo[1];
    holeFrac       = maskInfo[2];
    erodeUsed_um   = maskInfo[3];
    nFragments     = maskInfo[4];

    if (tissueArea_mm2 <= 0) {
        print("  FAILED: no tissue survived masking. Check TISSUE_OD_THRESHOLD.");
        return;
    }
    print("  tissue " + d2s(tissueArea_mm2,3) + " mm2 | holes " + d2s(holeFrac*100,1) +
          "% | frags " + nFragments + " | erode " + erodeUsed_um + " um");

    // ---------- 4. Colour deconvolution ----------
    selectWindow(src);
    run("Select None");
    run("Colour Deconvolution2", buildVectorString());

    c1 = findDeconWindow(src, 1);   // Hematoxylin
    c2 = findDeconWindow(src, 2);   // DAB / brown
    c3 = findDeconWindow(src, 3);   // Magenta
    if (c1 == "" || c2 == "" || c3 == "") {
        print("  FAILED: could not locate deconvolution output windows.");
        print("  Run Colour Deconvolution2 once from the GUI with the macro");
        print("  Recorder on, and compare its option string to buildVectorString().");
        return;
    }

    // ---------- 5a. SCAN MODE ----------
    if (THRESHOLD_SCAN) {
        if (SCAN_CHANNEL != "dab")
            emitThresholdScan(c3, "MASK_UNFILLED", base, outDir, tissueArea_mm2, "", "magenta_area_fraction", SCAN_MAX);
        if (SCAN_CHANNEL != "magenta")
            emitThresholdScan(c2, "MASK_UNFILLED", base, outDir, tissueArea_mm2, "dab_", "dab_area_fraction", DAB_SCAN_MAX);
        print("  scan done in " + d2s((getTime()-t0)/1000.0,1) + " s");
        return;
    }

    // ---------- 5b. PRIMARY: whole-section area fraction ----------
    magAF = areaFractionThroughMask(c3, "MASK_UNFILLED", MAGENTA_THRESHOLD);
    brnAF = areaFractionThroughMask(c2, "MASK_UNFILLED", BROWN_THRESHOLD);
    print("  magenta area fraction = " + d2s(magAF*100,3) + " %   <-- PRIMARY");
    print("  brown   area fraction = " + d2s(brnAF*100,3) + " %");
    // INTENSITY metric: mean magenta OD over POSITIVE pixels only (mag > threshold,
    // inside tissue). High-OD regime -> crosstalk-negligible, stays positive (unlike
    // mean-over-all-tissue). Normalize downstream by hema_mean for cellularity.
    magPosMean = positiveMeanThroughMask(c3, "MASK_UNFILLED", MAGENTA_THRESHOLD);
    print("  magenta positive-pixel mean OD = " + d2s(magPosMean,4));

    // ---------- 5c. SECONDARY: means (inspection only) ----------
    hU = measureThroughMask(c1, "MASK_UNFILLED");
    dU = measureThroughMask(c2, "MASK_UNFILLED");
    mU = measureThroughMask(c3, "MASK_UNFILLED");
    hF = measureThroughMask(c1, "MASK_FILLED");
    dF = measureThroughMask(c2, "MASK_FILLED");
    mF = measureThroughMask(c3, "MASK_FILLED");

    // ---------- 5d. Per-FOV grid ----------
    // FOV outputs (per-FOV table + numbered overlay) go in their own subfolder,
    // separate from the whole-sample outputs.
    fovDir = outDir + "fov" + File.separator;
    File.makeDirectory(fovDir);
    fovCsv = fovDir + base + "_fov.csv";
    nFov = emitFovGrid(c2, c3, "MASK_UNFILLED", W, H, base, fovCsv);
    print("  FOVs: " + nFov + " (grid " + FOV_GRID_UM + " um, FOV " + FOV_SIZE_PX +
          " px = " + d2s(FOV_SIZE_PX*PIXEL_UM,0) + " um)");
    if (nFov == 0) {
        print("  WARNING: no FOV fit entirely inside the eroded mask. Sample is small");
        print("  relative to FOV_SIZE_PX, or the grid pitch is too coarse. n=0 for this");
        print("  image -- it contributes to whole-section AF only.");
    }

    // ---------- 6. Per-image CSV ----------
    csv = outDir + base + "_absorbance.csv";
    f = File.open(csv);
    print(f, "image,pixel_um,I0_R,I0_G,I0_B,erode_um,n_fragments," +
             "tissue_area_mm2_unfilled,tissue_area_mm2_filled,hole_fraction," +
             "magenta_threshold,magenta_area_fraction,brown_threshold,brown_area_fraction," +
             "n_fov,fov_size_px,fov_grid_um," +
             "hema_mean_unfilled,hema_sd_unfilled,dab_mean_unfilled,dab_sd_unfilled," +
             "mag_mean_unfilled,mag_sd_unfilled," +
             "hema_mean_filled,hema_sd_filled,dab_mean_filled,dab_sd_filled," +
             "mag_mean_filled,mag_sd_filled," +
             "threshold_placeholder,vector_DdotM,mag_pos_mean");
    print(f, base + "," + PIXEL_UM + "," +
             d2s(i0[0],3) + "," + d2s(i0[1],3) + "," + d2s(i0[2],3) + "," +
             erodeUsed_um + "," + nFragments + "," +
             d2s(tissueArea_mm2,5) + "," + d2s(filledArea_mm2,5) + "," + d2s(holeFrac,5) + "," +
             d2s(MAGENTA_THRESHOLD,4) + "," + d2s(magAF,6) + "," +
             d2s(BROWN_THRESHOLD,4) + "," + d2s(brnAF,6) + "," +
             nFov + "," + FOV_SIZE_PX + "," + FOV_GRID_UM + "," +
             d2s(hU[0],6) + "," + d2s(hU[1],6) + "," + d2s(dU[0],6) + "," + d2s(dU[1],6) + "," +
             d2s(mU[0],6) + "," + d2s(mU[1],6) + "," +
             d2s(hF[0],6) + "," + d2s(hF[1],6) + "," + d2s(dF[0],6) + "," + d2s(dF[1],6) + "," +
             d2s(mF[0],6) + "," + d2s(mF[1],6) + "," +
             THRESHOLD_IS_PLACEHOLDER + "," + VECTOR_DdotM + "," + d2s(magPosMean,6));
    File.close(f);
    print("  CSV -> " + csv);

    // ---------- 7. QC + channel outputs ----------
    if (SAVE_THRESH_QC) saveThresholdQC(c3, W, H, base, outDir);

    if (SAVE_CHANNEL_TIFS || SAVE_FIGURE_PNGS) {
        writeChannel(c1, outDir + base + "_1_hematoxylin",    COL_H_R, COL_H_G, COL_H_B, LUT_H_MIN, LUT_H_MAX);
        writeChannel(c2, outDir + base + "_2_dab_IL15",       COL_D_R, COL_D_G, COL_D_B, LUT_D_MIN, LUT_D_MAX);
        writeChannel(c3, outDir + base + "_3_magenta_IL15Ra", COL_M_R, COL_M_G, COL_M_B, LUT_M_MIN, LUT_M_MAX);
    }

    if (SAVE_STAIN_FIGS) saveChannelFigures(src, c1, c2, c3, base, outDir);

    if (SAVE_FOV_OVERLAY) saveFovOverlay(src, base, fovDir, fovCsv, nFov, W);

    print("  done in " + d2s((getTime()-t0)/1000.0,1) + " s");
}


// ============================================================
// =================== THRESHOLD SCAN =========================
// ============================================================

function emitThresholdScan(chanWin, maskWin, base, outDir, tissueArea_mm2, tag, colName, scanMax) {
    // tag "" with colName "magenta_area_fraction" gives the Vector Red scan's file names.
    csv = outDir + base + "_" + tag + "threshold_scan.csv";
    f = File.open(csv);
    print(f, "image,threshold," + colName + ",positive_area_mm2");

    print("  threshold   area_fraction");
    for (t = SCAN_MIN; t <= scanMax + 1e-9; t += SCAN_STEP) {
        af = areaFractionThroughMask(chanWin, maskWin, t);
        print(f, base + "," + d2s(t,4) + "," + d2s(af,6) + "," + d2s(af*tissueArea_mm2,6));
        print("     " + d2s(t,2) + "        " + d2s(af*100,3) + " %");
    }
    File.close(f);
    print("  scan -> " + csv);

    if (SAVE_THRESH_QC) {
        for (t = SCAN_MIN; t <= scanMax + 1e-9; t += SCAN_STEP*2) {
            saveScanQC(chanWin, maskWin, base, outDir, t, tag);
        }
        print("  overlay PNGs -> " + outDir + base + "_" + tag + "scan_thr*.png");
    }
}

function saveScanQC(chanWin, maskWin, base, outDir, thr, tag) {
    selectWindow(chanWin);
    run("Select None");
    run("Duplicate...", "title=__scanq");
    setThreshold(thr, 1e30);
    run("Convert to Mask");
    imageCalculator("AND", "__scanq", maskWin);

    selectWindow("__scanq");
    getDimensions(w, h, c, s, fr);
    qw = 1400;
    if (w < qw) qw = w;
    qh = round(h * qw / w);
    run("Size...", "width=" + qw + " height=" + qh + " depth=1 interpolation=None");
    saveAs("PNG", outDir + base + "_" + tag + "scan_thr" + d2s(thr,2) + ".png");
    close();
}


// ============================================================
// ======================= FOV GRID ===========================
// ============================================================

function emitFovGrid(brownWin, magWin, maskWin, W, H, base, csvPath) {
    // Systematic uniform random sampling: single random start, fixed pitch.
    // Better variance properties than independent random placement, and it is
    // the stereology-supported approach.
    pitch = round(FOV_GRID_UM / PIXEL_UM);
    if (pitch < FOV_SIZE_PX) {
        print("  NOTE: grid pitch (" + pitch + " px) < FOV size (" + FOV_SIZE_PX +
              " px) -- FOVs overlap; samples are not independent.");
    }

    seed = FOV_SEED_BASE + strHash(base);
    random("seed", seed);
    x0 = floor(random * pitch);
    y0 = floor(random * pitch);

    f = File.open(csvPath);
    print(f, "image,fov_id,x_px,y_px,fov_size_px,tissue_frac_in_fov," +
             "magenta_area_fraction,brown_area_fraction,mag_mean,brown_mean,placement");

    n = 0;
    px = newArray(0); py = newArray(0);    // placed FOV top-left coords (for top-up overlap test)
    for (y = y0; y + FOV_SIZE_PX <= H; y += pitch) {
        for (x = x0; x + FOV_SIZE_PX <= W; x += pitch) {

            selectWindow(maskWin);
            makeRectangle(x, y, FOV_SIZE_PX, FOV_SIZE_PX);
            getStatistics(a, mMean);
            tfrac = mMean / 255.0;
            run("Select None");

            if (FOV_REQUIRE_FULL) {
                if (tfrac < FOV_MIN_TISSUE_FRAC) continue;
            } else {
                if (tfrac <= 0) continue;
            }

            magAF = areaFractionInRect(magWin,   maskWin, x, y, FOV_SIZE_PX, MAGENTA_THRESHOLD);
            brnAF = areaFractionInRect(brownWin, maskWin, x, y, FOV_SIZE_PX, BROWN_THRESHOLD);

            selectWindow(magWin);
            makeRectangle(x, y, FOV_SIZE_PX, FOV_SIZE_PX);
            getStatistics(a2, magMean);
            run("Select None");

            selectWindow(brownWin);
            makeRectangle(x, y, FOV_SIZE_PX, FOV_SIZE_PX);
            getStatistics(a3, brnMean);
            run("Select None");

            print(f, base + "," + n + "," + x + "," + y + "," + FOV_SIZE_PX + "," +
                     d2s(tfrac,5) + "," + d2s(magAF,6) + "," + d2s(brnAF,6) + "," +
                     d2s(magMean,6) + "," + d2s(brnMean,6) + ",grid");
            px = Array.concat(px, x); py = Array.concat(py, y);
            n++;
        }
    }
    nGrid = n;

    // ---- Minimum-FOV top-up ----
    // The systematic grid can under-sample thin/curved sections (a full FOV rarely
    // lands entirely inside the eroded strip at the fixed pitch). If it placed fewer
    // than MIN_FOV_PER_SAMPLE, scan the eroded mask on a non-overlapping step and add
    // full-size FOVs that lie ENTIRELY inside tissue and do not overlap an already-
    // placed FOV, until the target is met or none fit. Flagged placement=topup so the
    // systematic-grid FOVs remain identifiable as the primary unbiased sample.
    if (MIN_FOV_PER_SAMPLE > 0 && n < MIN_FOV_PER_SAMPLE) {
        step = maxOf(round(FOV_SIZE_PX/2), 1);   // half-FOV scan: finer placement on
                                                 // curved/thin strips; non-overlap still
                                                 // enforced below, so FOVs stay >=1 FOV apart
        for (yy = 0; yy + FOV_SIZE_PX <= H && n < MIN_FOV_PER_SAMPLE; yy += step) {
            for (xx = 0; xx + FOV_SIZE_PX <= W && n < MIN_FOV_PER_SAMPLE; xx += step) {
                selectWindow(maskWin);
                makeRectangle(xx, yy, FOV_SIZE_PX, FOV_SIZE_PX);
                getStatistics(aT, mT);
                run("Select None");
                if (mT/255.0 < 0.999) continue;                 // must be fully inside tissue
                ov = false;                                      // reject overlap with placed FOVs
                for (k = 0; k < px.length; k++)
                    if (abs(xx - px[k]) < FOV_SIZE_PX && abs(yy - py[k]) < FOV_SIZE_PX) ov = true;
                if (ov) continue;
                magAF = areaFractionInRect(magWin,   maskWin, xx, yy, FOV_SIZE_PX, MAGENTA_THRESHOLD);
                brnAF = areaFractionInRect(brownWin, maskWin, xx, yy, FOV_SIZE_PX, BROWN_THRESHOLD);
                selectWindow(magWin);   makeRectangle(xx,yy,FOV_SIZE_PX,FOV_SIZE_PX); getStatistics(a2,magMean); run("Select None");
                selectWindow(brownWin); makeRectangle(xx,yy,FOV_SIZE_PX,FOV_SIZE_PX); getStatistics(a3,brnMean); run("Select None");
                print(f, base + "," + n + "," + xx + "," + yy + "," + FOV_SIZE_PX + "," +
                         d2s(mT/255.0,5) + "," + d2s(magAF,6) + "," + d2s(brnAF,6) + "," +
                         d2s(magMean,6) + "," + d2s(brnMean,6) + ",topup");
                px = Array.concat(px, xx); py = Array.concat(py, yy);
                n++;
            }
        }
        if (n > nGrid)
            print("  FOV top-up: grid " + nGrid + " -> +" + (n - nGrid) +
                  " non-overlapping full FOVs (target " + MIN_FOV_PER_SAMPLE + ")");
    }

    // ---- Adaptive fallback: guarantee >= 1 FOV per sample ----
    // Small/thin sections that cannot fully contain any grid FOV would yield
    // n_fov = 0. Place ONE FOV at the deepest-interior point of the eroded mask
    // (max of its distance map), sized to the largest square that fits there
    // (capped at FOV_SIZE_PX). Its true size goes in fov_size_px so a fallback
    // FOV is never silently pooled with full-size FOVs downstream.
    if (n == 0) {
        selectWindow(maskWin);
        run("Select None");
        run("Duplicate...", "title=__fbm");
        dsf = 8;
        dwf = maxOf(round(W/dsf), 4);
        dhf = maxOf(round(H/dsf), 4);
        run("Size...", "width=" + dwf + " height=" + dhf + " depth=1 interpolation=None");
        setThreshold(128, 255);
        run("Convert to Mask");
        run("Distance Map");
        getStatistics(aFb, mFb, mnFb, mxFb);       // mxFb = max clearance (downsampled px)
        fbSize = round(2 * mxFb * dsf * 0.9);       // largest square that safely fits (full-res px)
        if (fbSize > FOV_SIZE_PX) fbSize = FOV_SIZE_PX;
        if (fbSize >= 64) {                          // skip true slivers (< ~48 um)
            setThreshold(mxFb - 0.5, 255);
            run("Create Selection");
            getSelectionBounds(bxf, byf, bwf, bhf);
            close("__fbm");
            fxF = round((bxf + bwf/2) * dsf) - floor(fbSize/2);
            fyF = round((byf + bhf/2) * dsf) - floor(fbSize/2);
            if (fxF < 0) fxF = 0;
            if (fyF < 0) fyF = 0;
            if (fxF + fbSize > W) fxF = W - fbSize;
            if (fyF + fbSize > H) fyF = H - fbSize;
            selectWindow(maskWin);
            makeRectangle(fxF, fyF, fbSize, fbSize);
            getStatistics(aT, tMeanF);
            tfracF = tMeanF / 255.0;
            run("Select None");
            magAFf = areaFractionInRect(magWin,   maskWin, fxF, fyF, fbSize, MAGENTA_THRESHOLD);
            brnAFf = areaFractionInRect(brownWin, maskWin, fxF, fyF, fbSize, BROWN_THRESHOLD);
            selectWindow(magWin);   makeRectangle(fxF, fyF, fbSize, fbSize); getStatistics(a2f, magMeanF); run("Select None");
            selectWindow(brownWin); makeRectangle(fxF, fyF, fbSize, fbSize); getStatistics(a3f, brnMeanF); run("Select None");
            print(f, base + ",0," + fxF + "," + fyF + "," + fbSize + "," +
                     d2s(tfracF,5) + "," + d2s(magAFf,6) + "," + d2s(brnAFf,6) + "," +
                     d2s(magMeanF,6) + "," + d2s(brnMeanF,6) + ",fallback");
            n = 1;
            print("  FOV fallback: 0 grid FOVs -> placed 1 (" + fbSize + " px = " +
                  d2s(fbSize*PIXEL_UM,0) + " um, tissue_frac " + d2s(tfracF,3) + ")");
        } else {
            close("__fbm");
            print("  FOV fallback: tissue too thin for any FOV; n_fov stays 0.");
        }
    }

    File.close(f);
    return n;
}

function areaFractionInRect(chanWin, maskWin, x, y, sz, thr) {
    // Positive pixels in the FOV as a fraction of TISSUE pixels in the FOV.
    // sz = FOV side in px (variable so the adaptive fallback can use a smaller box).
    selectWindow(chanWin);
    makeRectangle(x, y, sz, sz);
    run("Duplicate...", "title=__fov_c");
    run("Select None");
    setThreshold(thr, 1e30);
    run("Convert to Mask");

    selectWindow(maskWin);
    makeRectangle(x, y, sz, sz);
    run("Duplicate...", "title=__fov_m");
    run("Select None");
    setThreshold(128, 255);
    run("Convert to Mask");

    imageCalculator("AND", "__fov_c", "__fov_m");

    selectWindow("__fov_c");
    getStatistics(a, posMean);
    nPos = (posMean/255.0) * sz * sz;

    selectWindow("__fov_m");
    getStatistics(a2, tisMean);
    nTis = (tisMean/255.0) * sz * sz;

    close("__fov_c");
    close("__fov_m");
    selectWindow(chanWin); run("Select None");
    selectWindow(maskWin); run("Select None");

    if (nTis <= 0) return 0;
    return nPos / nTis;
}

function strHash(s) {
    h = 0;
    for (i = 0; i < lengthOf(s); i++) {
        h = (h * 31 + charCodeAt(s, i)) % 100000;
    }
    return h;
}


// ============================================================
// ====================== MEASUREMENT =========================
// ============================================================

function areaFractionThroughMask(chanWin, maskWin, thr) {
    // PRIMARY ENDPOINT. Fraction of TISSUE pixels above thr.
    // Keys on high-contrast discrete signal; avoids the diffuse low-OD regime
    // where the ill-conditioned inversion puts most of its error.
    selectWindow(chanWin);
    run("Select None");
    run("Duplicate...", "title=__af_c");
    setThreshold(thr, 1e30);
    run("Convert to Mask");

    imageCalculator("AND", "__af_c", maskWin);

    selectWindow("__af_c");
    run("Select None");
    getStatistics(a, posMean);
    getDimensions(w, h, c, s, fr);
    nPos = (posMean/255.0) * w * h;
    close("__af_c");

    selectWindow(maskWin);
    run("Select None");
    getStatistics(a2, tisMean);
    getDimensions(w2, h2, c2, s2, f2);
    nTis = (tisMean/255.0) * w2 * h2;

    if (nTis <= 0) return 0;
    return nPos / nTis;
}

function measureThroughMask(chanWin, maskWin) {
    // Inspection only. Not a valid endpoint with near-collinear vectors (see header).
    selectWindow(maskWin);
    run("Select None");
    setThreshold(128, 255);
    run("Create Selection");
    if (selectionType() == -1) { r = newArray(0, 0); return r; }
    roiManager("reset");
    roiManager("Add");

    selectWindow(chanWin);
    roiManager("Select", 0);
    getStatistics(area, mean, mn, mx, sd);
    run("Select None");
    roiManager("reset");

    r = newArray(mean, sd);
    return r;
}

function positiveMeanThroughMask(chanWin, maskWin, thr) {
    // INTENSITY of the positive signal: mean chan OD over pixels where chan > thr
    // AND inside the tissue mask. Uses only crash-safe ops (8-bit Convert-to-Mask,
    // 8-bit imageCalculator AND, Create Selection + getStatistics MEAN -- the same
    // pattern as measureThroughMask). Does NOT use getValue("Median") or an
    // imageCalculator multiply on 32-bit, because both hang under xvfb.
    selectWindow(chanWin);
    run("Select None");
    run("Duplicate...", "title=__pm");
    setThreshold(thr, 1e30);
    run("Convert to Mask");                 // 8-bit: 255 where chan > thr
    imageCalculator("AND", "__pm", maskWin); // positive AND tissue (8-bit AND -- safe)
    selectWindow("__pm");
    setThreshold(128, 255);
    run("Create Selection");
    if (selectionType() == -1) { close("__pm"); return 0; }   // no positive pixels
    roiManager("reset");
    roiManager("Add");
    selectWindow(chanWin);
    roiManager("Select", 0);
    getStatistics(area, mean, mn, mx, sd);
    run("Select None");
    roiManager("reset");
    close("__pm");
    return mean;
}


// ============================================================
// ==================== WHITE REFERENCE =======================
// ============================================================

function measureI0(src, W, H) {
    out = newArray(3);
    if (I0_MANUAL) { out[0]=I0_R; out[1]=I0_G; out[2]=I0_B; return out; }

    selectWindow(src);
    pw = minOf(round(W*0.04), 400);
    ph = minOf(round(H*0.04), 400);
    if (pw < 10) pw = 10;
    if (ph < 10) ph = 10;

    xs = newArray(0, W-pw, 0,    W-pw);
    ys = newArray(0, 0,    H-ph, H-ph);

    bestLum = -1;
    for (k = 0; k < 4; k++) {
        selectWindow(src);
        makeRectangle(xs[k], ys[k], pw, ph);
        rgb = meanRGBofSelection(src);
        lum = 0.299*rgb[0] + 0.587*rgb[1] + 0.114*rgb[2];
        if (lum > bestLum) { bestLum = lum; out[0]=rgb[0]; out[1]=rgb[1]; out[2]=rgb[2]; }
    }
    selectWindow(src);
    run("Select None");

    if (bestLum < 150) {
        print("  WARNING: brightest corner luminance is " + d2s(bestLum,1) +
              " -- corners may not be clean glass. Consider I0_MANUAL=1.");
    }
    for (k = 0; k < 3; k++) if (out[k] < 1) out[k] = 1;
    return out;
}

function meanRGBofSelection(src) {
    selectWindow(src);
    run("Duplicate...", "title=__i0_tmp");
    run("RGB Stack");
    o = newArray(3);
    for (c = 1; c <= 3; c++) { setSlice(c); getStatistics(a, mn); o[c-1] = mn; }
    close("__i0_tmp");
    return o;
}


// ============================================================
// ====================== TISSUE MASK =========================
// ============================================================

function buildMasks(src, W, H, i0, base, outDir) {
    ds = MASK_DOWNSAMPLE;
    dw = maxOf(round(W/ds), 4);
    dh = maxOf(round(H/ds), 4);
    dsPixelUm = PIXEL_UM * ds;    // 0.75488 * 4 = 3.02 um per mask pixel

    selectWindow(src);
    run("Select None");
    run("Duplicate...", "title=__ds");
    run("Size...", "width=" + dw + " height=" + dh + " depth=1 average interpolation=Bilinear");

    // total OD = sum of per-channel OD, referenced to measured glass.
    // Stain-agnostic "how much is here" map -- detects faint hematoxylin-only
    // tissue that a saturation threshold would miss.
    selectWindow("__ds");
    run("RGB Stack");
    run("32-bit");
    for (c = 1; c <= 3; c++) {
        setSlice(c);
        run("Divide...", "value=" + i0[c-1] + " slice");
        run("Min...", "value=0.0001 slice");         // keep log finite
        run("Log", "slice");                          // natural log
        run("Divide...", "value=2.302585 slice");     // -> log10
        run("Multiply...", "value=-1 slice");         // -> OD
    }
    // 8-bit copy of the OD stack for the chromaticity gate (OD 0-0.5 -> 0-255).
    if (TISSUE_CHROMA_LEVELS > 0) {
        run("Duplicate...", "title=__odstack8 duplicate");
        selectWindow("__odstack8"); setMinAndMax(0, 0.5); run("8-bit");
        selectWindow("__ds");
    }

    // Max single-channel OD map, for the near-white glass gate below. Computed BEFORE the
    // sum projection because Z Project consumes the stack selection.
    run("Z Project...", "projection=[Max Intensity]");
    rename("__odmax");
    selectWindow("__ds");
    run("Z Project...", "projection=[Sum Slices]");
    rename("__odsum");
    close("__ds");

    // ---- Threshold on a SMOOTHED OD map ----
    // Faint/sparse staining fires the threshold on scattered pixels; a light
    // Gaussian blur on the OD sum reconnects them so thin tissue (biopsy strips,
    // rings, multi-segment sections) forms connected regions instead of shattering
    // into sub-threshold fragments that the size filter would then delete.
    selectWindow("__odsum");
    if (MASK_SMOOTH_PX > 0) run("Gaussian Blur...", "sigma=" + MASK_SMOOTH_PX);
    setThreshold(TISSUE_OD_THRESHOLD, 1e30);
    run("Convert to Mask");
    if (is("Inverting LUT")) run("Invert LUT");     // tissue = 255 (so the AND below is valid)
    rename("__m0");
    // ---- CHROMATICITY gate: glass is grey; faint tissue keeps a stain hue ----
    // All 8-bit: 8-bit imageCalculator is safe, the 32-bit form hangs under xvfb.
    if (TISSUE_CHROMA_LEVELS > 0) {
        selectWindow("__odstack8");
        run("Z Project...", "projection=[Max Intensity]"); rename("__qmax");
        selectWindow("__odstack8");
        run("Z Project...", "projection=[Min Intensity]"); rename("__qmin");
        imageCalculator("Subtract create", "__qmax", "__qmin"); rename("__chroma");
        setThreshold(TISSUE_CHROMA_LEVELS, 255);
        run("Convert to Mask");
        if (is("Inverting LUT")) run("Invert LUT");
        rename("__mchroma");
        imageCalculator("AND", "__m0", "__mchroma");
        close("__mchroma");
        if (isOpen("__qmax"))     close("__qmax");
        if (isOpen("__qmin"))     close("__qmin");
        if (isOpen("__chroma"))   close("__chroma");
        if (isOpen("__odstack8")) close("__odstack8");
        print("  chromaticity gate applied at " + TISSUE_CHROMA_LEVELS + " levels");
    }

    // ---- Near-white GLASS rejection (see TISSUE_MAXOD_THRESHOLD) ----
    // Require real absorbance in at least one channel. 8-bit AND is used because
    // imageCalculator on 32-bit hangs under xvfb.
    if (TISSUE_MAXOD_THRESHOLD > 0) {
        selectWindow("__odmax");
        if (MASK_SMOOTH_PX > 0) run("Gaussian Blur...", "sigma=" + MASK_SMOOTH_PX);
        setThreshold(TISSUE_MAXOD_THRESHOLD, 1e30);
        run("Convert to Mask");
        if (is("Inverting LUT")) run("Invert LUT");
        rename("__mmax");
        imageCalculator("AND", "__m0", "__mmax");
        close("__mmax");
    }
    if (isOpen("__odmax")) close("__odmax");

    run("Options...", "iterations=1 count=1 black");
    run("Open");                                   // remove single-pixel speckle

    run("Set Scale...", "distance=1 known=" + dsPixelUm + " unit=micron");
    minA = MIN_FRAGMENT_MM2 * 1e6;                 // mm^2 -> um^2
    run("Analyze Particles...", "size=" + minA + "-Infinity pixel show=Masks clear");
    rename("__mfrag");
    if (is("Inverting LUT")) run("Invert LUT");
    close("__m0");
    run("Set Scale...", "distance=1 known=" + dsPixelUm + " unit=micron");

    nFrag = countFragments("__mfrag", dsPixelUm, minA);

    // __m1_filled: ALL interior holes filled (the "filled" comparison denominator).
    selectWindow("__mfrag");
    run("Duplicate...", "title=__m1_filled");
    run("Fill Holes");

    // __m1: PRIMARY tissue mask.
    if (EXCLUDE_LARGE_LUMENS) {
        // Fill interior holes SMALLER than HOLE_FILL_MAXMM2 but keep genuinely large
        // lumens (vessels, empty bays) EXCLUDED. Built as (all-filled) minus (large holes).
        holeCapUm2 = HOLE_FILL_MAXMM2 * 1e6;
        imageCalculator("Subtract create", "__m1_filled", "__mfrag");   // interior holes
        rename("__holes");
        run("Set Scale...", "distance=1 known=" + dsPixelUm + " unit=micron");
        run("Analyze Particles...", "size=" + holeCapUm2 + "-Infinity pixel show=Masks clear");
        rename("__biglumens");                          // holes large enough to keep excluded
        if (is("Inverting LUT")) run("Invert LUT");
        imageCalculator("Subtract create", "__m1_filled", "__biglumens");
        rename("__m1");
        if (is("Inverting LUT")) run("Invert LUT");
        close("__holes"); close("__biglumens"); close("__mfrag");
    } else {
        // Primary mask = the hole-filled tissue silhouette (no large-lumen carve-out).
        selectWindow("__m1_filled");
        run("Duplicate...", "title=__m1");
        close("__mfrag");
    }

    // Edge erosion, with PROGRESSIVE relaxation for thin/small tissue. Thin core
    // biopsies are almost entirely edge: a fixed 100 um erosion (from both sides)
    // can remove the whole strip, which would zero the tissue area and drop the
    // sample. Start at EDGE_ERODE_UM; if it removes more than
    // EDGE_ERODE_MAXFRAC of the tissue -- or removes ALL of it -- step the erosion
    // down the ladder until enough tissue survives. The value actually used is
    // returned in the erode_um column, so erosion is transparent per image.
    areaBefore = maskAreaMM2("__m1", dsPixelUm);
    erodeLadder = newArray(EDGE_ERODE_UM, EDGE_ERODE_FALLBACK_UM,
                           round(EDGE_ERODE_FALLBACK_UM/2), 0);
    erodeUm = erodeLadder[0];
    erodedName = erodeMask("__m1", erodeUm, dsPixelUm, "__m1_er");
    areaAfter  = maskAreaMM2(erodedName, dsPixelUm);

    li = 1;
    while (areaBefore > 0 && li < erodeLadder.length &&
           (areaAfter <= 0 || (1.0 - areaAfter/areaBefore) > EDGE_ERODE_MAXFRAC)) {
        close(erodedName);
        erodeUm = erodeLadder[li];
        erodedName = erodeMask("__m1", erodeUm, dsPixelUm, "__m1_er");
        areaAfter  = maskAreaMM2(erodedName, dsPixelUm);
        li++;
    }
    if (erodeUm < EDGE_ERODE_UM) {
        print("  NOTE: thin/small tissue -- edge erosion relaxed from " + EDGE_ERODE_UM +
              " um to " + erodeUm + " um so tissue survives (recorded in erode_um column;");
        print("  erosion is therefore NOT constant across the cohort for such images).");
    }
    erodedFilled = erodeMask("__m1_filled", erodeUm, dsPixelUm, "__m1f_er");
    areaAfterF = maskAreaMM2(erodedFilled, dsPixelUm);

    holeFrac = 0;
    if (areaAfterF > 0) holeFrac = (areaAfterF - areaAfter) / areaAfterF;


    upscaleMask(erodedName,   W, H, "MASK_UNFILLED");
    upscaleMask(erodedFilled, W, H, "MASK_FILLED");
    upscaleMask("__m1",       W, H, "MASK_TISSUE");   // pre-erosion cellular tissue (QC)

    if (SAVE_MASK_QC) saveMaskQC(src, W, H, base, outDir);

    close("__m1"); close("__m1_filled");
    close(erodedName); close(erodedFilled);
    if (isOpen("MASK_TISSUE")) close("MASK_TISSUE");
    if (isOpen("__odsum")) close("__odsum");

    r = newArray(areaAfter, areaAfterF, holeFrac, erodeUm, nFrag);
    return r;
}

function erodeMask(maskWin, um, dsPixelUm, newName) {
    selectWindow(maskWin);
    run("Duplicate...", "title=" + newName);
    px = round(um / dsPixelUm);    // 100 / 3.02 = 33 mask px -> 99.6 um actual
    if (px < 1) return newName;
    // Distance-map erosion: exact, and far faster than 33 Erode iterations.
    run("Duplicate...", "title=__edt_tmp");
    run("Distance Map");
    setThreshold(px, 1e30);
    run("Convert to Mask");
    rename("__edt_res");
    close(newName);
    selectWindow("__edt_res");
    rename(newName);
    if (isOpen("__edt_tmp")) close("__edt_tmp");
    return newName;
}

function upscaleMask(maskWin, W, H, outName) {
    selectWindow(maskWin);
    run("Duplicate...", "title=" + outName);
    // interpolation=None keeps the mask strictly binary on upscale.
    run("Size...", "width=" + W + " height=" + H + " depth=1 interpolation=None");
    setThreshold(128, 255);
    run("Convert to Mask");
    run("Set Scale...", "distance=1 known=" + PIXEL_UM + " unit=micron");
}

function maskAreaMM2(maskWin, dsPixelUm) {
    selectWindow(maskWin);
    run("Select None");
    getStatistics(area, mean);
    getDimensions(w, h, c, s, f);
    nWhite = (mean/255.0) * w * h;
    return nWhite * dsPixelUm * dsPixelUm / 1e6;
}

function countFragments(maskWin, dsPixelUm, minA_um2) {
    selectWindow(maskWin);
    run("Select None");
    run("Set Scale...", "distance=1 known=" + dsPixelUm + " unit=micron");
    run("Set Measurements...", "area redirect=None decimal=4");
    run("Analyze Particles...", "size=" + minA_um2 + "-Infinity pixel clear");
    n = nResults;
    run("Clear Results");
    return n;
}

function saveMaskQC(src, W, H, base, outDir) {
    // Downsampled RGB with TWO boundaries drawn on it:
    //   cyan   = tissue boundary before erosion
    //   yellow = the eroded analysis mask -- exactly what every area fraction is
    //            measured through.
    // The cyan->yellow band is the ~100 um edge rim that erosion removes. LOOK AT
    // THESE for the frayed / thin samples before trusting any number.
    qw = 1400;
    if (W < qw) qw = W;
    qh = round(H * qw / W);

    selectWindow(src);
    run("Select None");
    run("Duplicate...", "title=__qc");
    run("Size...", "width=" + qw + " height=" + qh + " depth=1 average interpolation=Bilinear");

    if (isOpen("MASK_TISSUE"))
        addMaskOutline("MASK_TISSUE", "__qc", qw, qh, "cyan");
    addMaskOutline("MASK_UNFILLED", "__qc", qw, qh, "yellow");

    selectWindow("__qc");
    run("Select None");
    run("Flatten");
    saveAs("PNG", outDir + base + "_maskQC.png");
    close();
    if (isOpen("__qc")) close("__qc");
}

// Draw the boundary of a binary mask onto the QC image as a coloured overlay.
function addMaskOutline(maskWin, qcWin, qw, qh, colorName) {
    selectWindow(maskWin);
    run("Select None");
    run("Duplicate...", "title=__qcmask");
    run("Size...", "width=" + qw + " height=" + qh + " depth=1 interpolation=None");
    setThreshold(128, 255);
    run("Create Selection");
    if (selectionType() != -1) {
        roiManager("reset");
        roiManager("Add");
        selectWindow(qcWin);
        roiManager("Select", 0);
        run("Properties... ", "stroke=" + colorName + " width=3");
        run("Add Selection...");
        roiManager("reset");
    }
    close("__qcmask");
}

function saveThresholdQC(magWin, W, H, base, outDir) {
    // Magenta-positive pixels at the frozen threshold, restricted to tissue.
    // Compare against the source. This is how you catch a threshold that is
    // picking up brown crosstalk rather than real magenta.
    selectWindow(magWin);
    run("Select None");
    run("Duplicate...", "title=__tq");
    setThreshold(MAGENTA_THRESHOLD, 1e30);
    run("Convert to Mask");
    imageCalculator("AND", "__tq", "MASK_UNFILLED");

    selectWindow("__tq");
    qw = 1400;
    if (W < qw) qw = W;
    qh = round(H * qw / W);
    run("Size...", "width=" + qw + " height=" + qh + " depth=1 interpolation=None");
    saveAs("PNG", outDir + base + "_magentaPositive.png");
    close();
}


// ============================================================
// =================== DECONVOLUTION ==========================
// ============================================================

function buildVectorString() {
    // Colour_1 = Hematoxylin, Colour_2 = DAB/brown, Colour_3 = Magenta
    // If this errors: run Colour Deconvolution2 once from the GUI with the
    // macro Recorder on and copy the exact option string it captures. The
    // parameter keys have drifted across plugin releases.
    s = "vectors=[User values] output=[32bit_Absorbance] simulated hide" +
        " [r1]=" + d2s(VEC_H_R,5) + " [g1]=" + d2s(VEC_H_G,5) + " [b1]=" + d2s(VEC_H_B,5) +
        " [r2]=" + d2s(VEC_D_R,5) + " [g2]=" + d2s(VEC_D_G,5) + " [b2]=" + d2s(VEC_D_B,5) +
        " [r3]=" + d2s(VEC_M_R,5) + " [g3]=" + d2s(VEC_M_G,5) + " [b3]=" + d2s(VEC_M_B,5);
    return s;
}

function findDeconWindow(src, idx) {
    cands = newArray(
        src + "-(Colour_" + idx + ")A",
        src + " (RGB)-(Colour_" + idx + ")A",
        src + "-(Colour_" + idx + ")",
        src + " (RGB)-(Colour_" + idx + ")"
    );
    for (k = 0; k < cands.length; k++) if (isOpen(cands[k])) return cands[k];
    for (i = 1; i <= nImages; i++) {
        selectImage(i);
        t = getTitle();
        if (indexOf(t, "Colour_" + idx) >= 0) return t;
    }
    return "";
}


// ============================================================
// ======================== OUTPUT ============================
// ============================================================

function writeChannel(winName, outPathNoExt, rEnd, gEnd, bEnd, lutMin, lutMax) {
    selectWindow(winName);
    run("Select None");
    // White at zero absorbance -> stain colour at high absorbance.
    // Display transform only; the 32-bit values are untouched.
    reds = newArray(256); greens = newArray(256); blues = newArray(256);
    for (i = 0; i < 256; i++) {
        f = i / 255.0;
        reds[i]   = round(255*(1-f) + rEnd*f);
        greens[i] = round(255*(1-f) + gEnd*f);
        blues[i]  = round(255*(1-f) + bEnd*f);
    }
    setMinAndMax(lutMin, lutMax);
    setLut(reds, greens, blues);

    if (SAVE_CHANNEL_TIFS) saveAs("Tiff", outPathNoExt + ".tif");
    if (SAVE_FIGURE_PNGS) {
        run("Flatten");    // GUI only -- needs a canvas
        saveAs("PNG", outPathNoExt + ".png");
        close();
    }
}

// ---- Stain figures: Original RGB + colour PNG per channel + labelled montage ----
// Stain panels bake a white->stain-colour gradient into 24-bit RGB (via RGB
// Color, NOT Flatten) and are then MASKED TO TISSUE: off-tissue pixels are set
// to white so the figure shows exactly what is measured (glass tint, which is
// excluded from all reported numbers, is not shown). The Original panel is the
// raw RGB source for direct comparison and is intentionally NOT masked.

function makeSourcePanelRGB(srcWin, pw, outName) {
    selectWindow(srcWin);
    run("Select None");
    run("Duplicate...", "title=" + outName);
    getDimensions(w, h, c, s, fr);
    ph = round(h * pw / w);
    run("Size...", "width=" + pw + " height=" + ph + " depth=1 interpolation=Bilinear");
    if (bitDepth() != 24) run("RGB Color");
    return outName;
}

function makeStainPanelRGB(chanWin, maskWin, pw, lutMin, lutMax, rEnd, gEnd, bEnd, outName) {
    selectWindow(chanWin);
    run("Select None");
    run("Duplicate...", "title=" + outName);
    getDimensions(w, h, c, s, fr);
    ph = round(h * pw / w);
    run("Size...", "width=" + pw + " height=" + ph + " depth=1 interpolation=Bilinear");
    reds = newArray(256); greens = newArray(256); blues = newArray(256);
    for (i = 0; i < 256; i++) {
        f = i / 255.0;
        reds[i]   = round(255*(1-f) + rEnd*f);
        greens[i] = round(255*(1-f) + gEnd*f);
        blues[i]  = round(255*(1-f) + bEnd*f);
    }
    setMinAndMax(lutMin, lutMax);
    setLut(reds, greens, blues);
    run("RGB Color");            // bake LUT + range into 24-bit RGB (in place)
    rename(outName);
    whiteOutBackground(outName, maskWin, pw, ph);
    return outName;
}

// Set every non-tissue pixel of an RGB panel to white, using the tissue mask.
function whiteOutBackground(panelName, maskWin, pw, ph) {
    selectWindow(maskWin);
    run("Select None");
    run("Duplicate...", "title=__pm");
    run("Size...", "width=" + pw + " height=" + ph + " depth=1 interpolation=None");
    setThreshold(128, 255);
    run("Convert to Mask");        // tissue = 255
    setThreshold(128, 255);
    run("Create Selection");       // selection = tissue
    if (selectionType() == -1) { close("__pm"); return; }
    roiManager("reset");
    roiManager("add");
    close("__pm");
    selectWindow(panelName);
    roiManager("select", 0);
    run("Make Inverse");           // -> non-tissue background
    setColor(255, 255, 255);
    fill();
    run("Select None");
    roiManager("reset");
}

function pastePanel(win, x, y) {
    selectWindow(win);
    getDimensions(pw, ph, c, s, fr);
    run("Select All"); run("Copy");
    selectWindow("__stainmont");
    makeRectangle(x, y, pw, ph);   // clipboard is pw x ph -> paste lands at (x,y)
    run("Paste");
    run("Select None");
}

function saveChannelFigures(srcWin, c1, c2, c3, base, outDir) {
    pw = STAIN_FIG_WIDTH;
    p0 = makeSourcePanelRGB(srcWin, pw, "__pan0");
    p1 = makeStainPanelRGB(c1, "MASK_UNFILLED", pw, LUT_H_MIN, LUT_H_MAX, COL_H_R, COL_H_G, COL_H_B, "__panH");
    p2 = makeStainPanelRGB(c2, "MASK_UNFILLED", pw, LUT_D_MIN, LUT_D_MAX, COL_D_R, COL_D_G, COL_D_B, "__panD");
    p3 = makeStainPanelRGB(c3, "MASK_UNFILLED", pw, LUT_M_MIN, LUT_M_MAX, COL_M_R, COL_M_G, COL_M_B, "__panM");

    // individual colour PNGs (rename back afterwards -- saveAs renames the window)
    selectWindow(p0); saveAs("PNG", outDir + base + "_ch0_original");        rename(p0);
    selectWindow(p1); saveAs("PNG", outDir + base + "_ch1_hematoxylin");     rename(p1);
    selectWindow(p2); saveAs("PNG", outDir + base + "_ch2_IL15_DAB");        rename(p2);
    selectWindow(p3); saveAs("PNG", outDir + base + "_ch3_IL15Ra_magenta");  rename(p3);

    // labelled 4-panel montage: Original | Hematoxylin | IL-15 | IL-15Ra
    selectWindow(p0); getDimensions(pwid, pht, c, s, fr);
    gap = 14; lab = 46;
    mw = 4*pwid + 5*gap;
    mh = pht + lab + gap;
    newImage("__stainmont", "RGB white", mw, mh, 1);
    pastePanel(p0, gap,               lab);
    pastePanel(p1, 2*gap + pwid,      lab);
    pastePanel(p2, 3*gap + 2*pwid,    lab);
    pastePanel(p3, 4*gap + 3*pwid,    lab);
    selectWindow("__stainmont");
    setFont("SansSerif", 26, "antialiased");
    setColor(25, 25, 25);
    drawString("Original (RGB)",     gap + 6,            lab - 14);
    drawString("Hematoxylin",        2*gap + pwid + 6,   lab - 14);
    drawString("IL-15 (DAB)",        3*gap + 2*pwid + 6, lab - 14);
    drawString("IL-15Ra (magenta)",  4*gap + 3*pwid + 6, lab - 14);
    saveAs("PNG", outDir + base + "_channels_montage");
    close();
    if (isOpen(p0)) close(p0);
    if (isOpen(p1)) close(p1);
    if (isOpen(p2)) close(p2);
    if (isOpen(p3)) close(p3);
}


// Whole-slide thumbnail with every FOV drawn as a numbered box. The number is
// the fov_id column in <base>_fov.csv, so each box on the slide correlates 1:1
// with a row of per-FOV measurements. Reads the CSV back (single source of
// truth) rather than re-deriving positions.
function saveFovOverlay(srcWin, base, outDir, csvPath, nFov, W) {
    pw = FOV_OVERLAY_WIDTH;
    p  = makeSourcePanelRGB(srcWin, pw, "__fovsrc");   // downsampled RGB slide
    getDimensions(pwid, pht, c, s, fr);
    scale = pwid / W;                                  // full-res px -> thumbnail px

    // draw boxes + numbers onto the thumbnail
    selectWindow("__fovsrc");
    setLineWidth(3);
    setFont("SansSerif", 13, "antialiased");   // must fit the 15 px label chip
    if (File.exists(csvPath)) {
        lines = split(File.openAsString(csvPath), "\n");
        for (i = 1; i < lines.length; i++) {           // row 0 is the header
            if (lengthOf(lines[i]) < 5) continue;
            col = split(lines[i], ",");
            if (col.length < 5) continue;
            fid = col[1];
            rx  = round(parseInt(col[2]) * scale);
            ry  = round(parseInt(col[3]) * scale);
            rsz = round(parseInt(col[4]) * scale);
            setColor(255, 238, 0);                     // yellow box outline
            drawRect(rx, ry, rsz, rsz);
            // The label sits outside the box, above its top-left corner, so the field
            // stays fully visible.
            lw = 6 + lengthOf(fid) * 8;                // compact chip
            lh = 15;
            ly = ry - lh - 1;                          // above the box
            if (ly < 0) ly = ry + rsz + 1;             // top row: drop below instead
            if (ly + lh > pht) ly = ry - lh - 1;
            setColor(255, 238, 0);  fillRect(rx, ly, lw, lh);
            setColor(0, 0, 0);      drawString(fid, rx + 3, ly + lh - 3);
        }
    }

    // add a caption band above the slide
    band = 44;
    newImage("__fovcanvas", "RGB white", pwid, pht + band, 1);
    selectWindow("__fovsrc"); run("Select All"); run("Copy");
    selectWindow("__fovcanvas"); makeRectangle(0, band, pwid, pht); run("Paste"); run("Select None");
    setFont("SansSerif", 24, "antialiased");
    setColor(25, 25, 25);
    drawString(base + "   -   " + nFov + " FOV(s), numbered per " + base + "_fov.csv",
               8, band - 14);
    saveAs("PNG", outDir + base + "_fov_overlay");
    close();
    if (isOpen("__fovsrc")) close("__fovsrc");
    setLineWidth(1);
}


// ============================================================
// ======================== UTILS =============================
// ============================================================

function isImageFile(name) {
    n = toLowerCase(name);
    return endsWith(n, ".tif")  || endsWith(n, ".tiff") ||
           endsWith(n, ".png")  || endsWith(n, ".jpg")  || endsWith(n, ".jpeg") ||
           endsWith(n, ".czi")  || endsWith(n, ".nd2")  || endsWith(n, ".lif");
}
