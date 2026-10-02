#!/bin/bash
#SBATCH --job-name=il15_ihc
#SBATCH --account=dalawson_lab           # UCI HPC3 lab account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G                         # ~0.7-3.8 GB RGB BigTIF -> tens of GB through 32-bit. 64 = headroom.
#SBATCH --time=02:00:00                   # per image; cohort mode ~2 min observed on 0.8-1.9 GB inputs
#SBATCH --output=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/il15_%A_%a.out
#SBATCH --error=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/il15_%A_%a.err
#SBATCH --array=1-74%20

# ============================================================
# IL-15 IHC pipeline -- one Slurm array task per image.
#
# Submit (full cohort):
#   ROOT=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis
#   N=$(wc -l < ${ROOT}/imagelist.txt)
#   sbatch --array=1-${N}%20 ${ROOT}/hpc/run_il15_array.sh
# Re-run only specific tasks (e.g. those missing a CSV):
#   sbatch --array=9,11,17,18,20,32,33,41,60,68,69 ${ROOT}/hpc/run_il15_array.sh
#
# TWO HARD REQUIREMENTS:
#   1. Colour Deconvolution2 builds an AWT dialog, so it CANNOT run under
#      --headless (it throws HeadlessException and writes nothing). Run it
#      under xvfb-run (virtual display) with NO --headless flag.
#   2. Keyence exports are TILED BigTIFFs, and ImageJ open() hangs on them.
#      The macro therefore opens them through Bio-Formats non-interactively.
#      Keep that path.
#
# The Fiji/xvfb launcher can exit nonzero even on success. Success is therefore
# judged by a per-image output CSV written during this job, and the launcher exit
# code is ignored.
# ============================================================

set -uo pipefail

# ---- CONFIGURE ----
# Paths come from config/site.env. Fiji lives on /pub scratch; after a reinstall, copy
# macros/colour_deconvolution2.jar into Fiji.app/plugins/.
ROOT="${IL15_ROOT:-/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis}"
source "${ROOT}/config/site.env"
# Fiji starts with HOME at the project root, as in every production run.
HOME="${ROOT}"
MACRO="${ROOT}/macros/il15_ihc_pipeline.ijm"
IMAGELIST="${ROOT}/imagelist.txt"
SUBDIRS="${ROOT}/subdirs.txt"             # line N = per-sample subfolder <ID>_<type> for image N
OUTDIR="${ROOT}/per_sample_out"
JAVA_MEM="$(( ${SLURM_MEM_PER_NODE:-64000} / 1000 - 6 ))g"   # track --mem alloc, leave ~6G OS headroom

mkdir -p "${ROOT}/logs" "${OUTDIR}"

if [[ ! -x "${FIJI}" ]]; then
    echo "ERROR: Fiji launcher not found/executable: ${FIJI}" >&2; exit 1
fi
if [[ ! -f "${IMAGELIST}" ]]; then
    echo "ERROR: ${IMAGELIST} not found." >&2; exit 1
fi

IMG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${IMAGELIST}")
if [[ -z "${IMG}" || ! -f "${IMG}" ]]; then
    echo "ERROR: no image at line ${SLURM_ARRAY_TASK_ID} of ${IMAGELIST}: '${IMG}'" >&2; exit 1
fi

BASE=$(basename "${IMG}"); BASE="${BASE%.*}"

# Per-sample output subfolder (<ID>_<tissue_type>), one folder per sample.
SUB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SUBDIRS}")
if [[ -z "${SUB}" ]]; then
    echo "ERROR: no subdir at line ${SLURM_ARRAY_TASK_ID} of ${SUBDIRS}" >&2; exit 1
fi
SAMPLE_OUT="${OUTDIR}/${SUB}"
mkdir -p "${SAMPLE_OUT}"
EXPECTED_CSV="${SAMPLE_OUT}/${BASE}_absorbance.csv"

echo "=========================================="
echo "Task ${SLURM_ARRAY_TASK_ID} | ${IMG} ($(du -h "${IMG}" | cut -f1)) | -> ${SUB} | node $(hostname) | $(date)"
echo "=========================================="

# xvfb-run provides the virtual display the plugin needs, so the launch omits --headless.
# -port0 disables ImageJ's single-instance handoff, so co-located array tasks, or a stale
# /tmp stub from a previous job on the node, do not collide. Without it Fiji throws
# java.rmi.ConnectException and the macro never runs. The rm below clears stale stubs.
JOB_START_EPOCH=$(date +%s)
rm -f /tmp/ImageJ-${USER}-*.stub 2>/dev/null || true
# Per-sample sampling overrides (grid pitch / erosion / min tissue fraction / FOV size).
# Small biopsy cores need a tighter grid than whole sections. Thin liver-met cores also
# take a smaller FOV. A blank FOV size keeps the cohort default. Table:
#   data/decisions/fov_sampling_overrides.csv  (sample_id,subfolder,area_mm2,fov_grid_um,
#                                               edge_erode_um,fov_min_tissue_frac,
#                                               fov_size_um,reason)
OVR="${ROOT}/data/decisions/fov_sampling_overrides.csv"
EXTRA=""
if [[ -f "${OVR}" ]]; then
    # cols: 1=sample_id 2=subfolder 3=area 4=grid_um 5=erode_um 6=min_tissue_frac 7=fov_size_um
    LINE=$(awk -F, -v s="${SUB}" 'NR>1 && $2==s {print $4","$5","$6","$7; exit}' "${OVR}")
    if [[ -n "${LINE}" ]]; then
        EXTRA="|$(echo "${LINE}"|cut -d, -f1)|$(echo "${LINE}"|cut -d, -f2)|$(echo "${LINE}"|cut -d, -f3)|$(echo "${LINE}"|cut -d, -f4)"
        echo "Sampling override for ${SUB}: pitch/erode/minfrac/fovsize = ${LINE}"
    fi
fi

# Disable ImageJ's single-instance handoff. With concurrent array tasks on one node, a new
# Fiji otherwise blocks in SingleInstance.sendArguments() and produces NO output.
# Each concurrent Fiji MUST get its own java.util.prefs root. Concurrent tasks otherwise
# contend on the shared lock in $HOME/.java/.userPrefs, and a losing task dies with
#   java.util.prefs.BackingStoreException: Couldn't get file lock
# even after it has done real work.
IJPRIV="${TMPDIR:-/tmp}/ij_${SLURM_ARRAY_TASK_ID:-0}_$$"
mkdir -p "$IJPRIV/prefs"
export _JAVA_OPTIONS="-Dnet.imagej.legacy.SingleInstance=false -Djava.io.tmpdir=$IJPRIV -Djava.util.prefs.userRoot=$IJPRIV/prefs -Djava.util.prefs.systemRoot=$IJPRIV/prefs"
mkdir -p "${TMPDIR:-/tmp}/ij_${SLURM_ARRAY_TASK_ID:-0}_$$"
xvfb-run -a "${FIJI}" --mem="${JAVA_MEM}" --console -port0 -macro "${MACRO}" "${IMG}|${SAMPLE_OUT}${EXTRA}"
LAUNCH_RC=$?

# Require the CSV to be NEWER than this job's start. A stale CSV from a previous run would
# otherwise mask a crash, and the task would report OK.
if [[ -f "${EXPECTED_CSV}" && -n "$(find "${EXPECTED_CSV}" -newermt "@${JOB_START_EPOCH}" 2>/dev/null)" ]]; then
    echo "OK: ${EXPECTED_CSV} written (launcher rc=${LAUNCH_RC}, ignored). $(date)"
    exit 0
else
    echo "FAIL: no ${EXPECTED_CSV} produced (launcher rc=${LAUNCH_RC}). $(date)" >&2
    exit 1
fi
