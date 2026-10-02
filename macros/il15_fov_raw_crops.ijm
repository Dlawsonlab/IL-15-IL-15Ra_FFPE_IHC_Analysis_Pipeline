// il15_fov_raw_crops.ijm
// Export the ORIGINAL RGB pixels of each requested FOV at NATIVE resolution, unlabelled.
//
// On the downsampled overview PNG a 662 px field is only ~17-20 px across, too coarse to
// tell a tear from a vessel lumen or a sinusoid. fov_quality_score.py and score_all_crops.py
// score these native-resolution crops instead.
//
// Arg: "imgPath|outDir|fovSpec"  where fovSpec = "id:x:y:sz;id:x:y:sz;..."
//   (x,y = top-left in FULL-RES px; sz = FOV size px; taken from <base>_fov.csv)
// Output: <outDir>/fov/raw/<base>_fovNN_raw.png   (NN = FOV id, 2-digit, native size)

setBatchMode(true);
arg = getArgument();
p = split(arg, "|");
imgPath = p[0];  outDir = p[1];  fovSpec = p[2];
if (!endsWith(outDir, File.separator)) outDir = outDir + File.separator;
rawDir = outDir + "fov" + File.separator + "raw" + File.separator;
File.makeDirectory(outDir + "fov");
File.makeDirectory(rawDir);
base = File.getNameWithoutExtension(imgPath);

// open exactly as the pipeline does -- group_files=false is REQUIRED for Keyence BigTIFFs
run("Bio-Formats Importer", "open=[" + imgPath + "] autoscale color_mode=Default " +
    "view=Hyperstack stack_order=XYCZT use_virtual_stack=false group_files=false");
if (nImages == 0) exit("failed to open " + imgPath);
src = getTitle();
if (bitDepth() != 24) {
    run("Stack to RGB");
    if (isOpen(src)) { selectWindow(src); close(); }
    src = getTitle();
}
selectWindow(src);
getDimensions(W, H, cc, ss, ff);

entries = split(fovSpec, ";");
n = 0;
for (i = 0; i < entries.length; i++) {
    e = String.trim(entries[i]);
    if (e == "") continue;
    f = split(e, ":");
    if (f.length < 4) continue;
    id = parseInt(f[0]);  x = parseInt(f[1]);  y = parseInt(f[2]);  sz = parseInt(f[3]);
    // clamp to the canvas so an edge FOV cannot abort the run
    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x + sz > W) sz = W - x;
    if (y + sz > H) sz = H - y;
    if (sz < 4) continue;

    selectWindow(src);
    makeRectangle(x, y, sz, sz);
    run("Duplicate...", "title=__crop");
    saveAs("PNG", rawDir + base + "_fov" + IJ.pad(id, 2) + "_raw.png");
    close();
    n = n + 1;
}
selectWindow(src); close();
print("wrote " + n + " native-resolution crops to " + rawDir);
setBatchMode(false);
eval("script", "System.exit(0);");
