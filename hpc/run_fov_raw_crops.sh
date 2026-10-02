#!/bin/bash
#SBATCH --job-name=il15_rawcrop
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/rawcrop_%A_%a.out
#SBATCH --error=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/rawcrop_%A_%a.err
# Export EVERY FOV of a sample as a native-resolution RGB crop, so the artefact metrics are
# computed at a resolution where they mean something (see il15_fov_raw_crops.ijm).
set -uo pipefail
ROOT="${IL15_ROOT:-/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis}"
source "${ROOT}/config/site.env"
# Fiji starts with HOME at the project root, as in every production run.
HOME="${ROOT}"
MACRO=${ROOT}/macros/il15_fov_raw_crops.ijm
IMG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${ROOT}/imagelist.txt)
SUB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${ROOT}/subdirs.txt)
OUTDIR=${ROOT}/per_sample_out/${SUB}
BASE=$(basename "$IMG"); BASE=${BASE%.*}
FOVCSV=${OUTDIR}/fov/${BASE}_fov.csv

echo "task ${SLURM_ARRAY_TASK_ID} | ${SUB} | $(date)"
if [[ ! -f "$FOVCSV" ]]; then echo "NO fov.csv at $FOVCSV" >&2; exit 1; fi

# Give each task its own Java prefs. Concurrent Fiji instances otherwise contend on the
# prefs file lock (BackingStoreException), which kills tasks at random.
IJPRIV="${TMPDIR:-/tmp}/ij_${SLURM_ARRAY_TASK_ID:-0}_$$"
mkdir -p "$IJPRIV/prefs"
export _JAVA_OPTIONS="-Dnet.imagej.legacy.SingleInstance=false -Djava.io.tmpdir=$IJPRIV -Djava.util.prefs.userRoot=$IJPRIV/prefs -Djava.util.prefs.systemRoot=$IJPRIV/prefs"

# Crop every field. Scoring the whole grid makes the flag rate comparable between samples
# and tissues.
SPEC=$(python3 - "$FOVCSV" <<'PY'
import csv, sys
out = []
for r in csv.DictReader(open(sys.argv[1])):
    out.append(f'{r["fov_id"]}:{r["x_px"]}:{r["y_px"]}:{r["fov_size_px"]}')
print(";".join(out))
PY
)
if [[ -z "$SPEC" ]]; then echo "empty spec" >&2; exit 1; fi
echo "fields to crop: $(echo "$SPEC" | tr ';' '\n' | grep -c .)"

rm -f /tmp/ImageJ-${USER}-*.stub 2>/dev/null || true
# Launch with -batch. It prints macro errors to stdout. Under xvfb, -macro buries the same
# errors in an invisible modal dialog, and the job hangs at 0% CPU until the time limit.
xvfb-run -a "$FIJI" --mem=26g --console -port0 -batch "$MACRO" "${IMG}|${OUTDIR}|${SPEC}"
echo "rc=$? $(date)"
echo "=== crops written ==="
ls -1 "${OUTDIR}/fov/raw/"*_raw.png 2>/dev/null | wc -l
