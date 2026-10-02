#!/bin/bash
#SBATCH --job-name=il15_dabcheck
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:50:00
#SBATCH --output=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/dabcheck_%A_%a.out
#SBATCH --error=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/dabcheck_%A_%a.err
#
# Check the IL-15 (DAB) positivity threshold, BROWN_THRESHOLD = 0.30, against the stain
# controls, and scan DAB area across thresholds for the group comparison.
#
# Input lines are "label|image[|I0]". The default list, inputs_lists/dab_check_inputs.txt,
# holds the 22 sections of the Vector Red scan and then the controls: secondary-only brain
# and liver (no primary antibody) and the haematoxylin-, DAB- and Vector Red-only images.
# The single-stain images have no glass at the corners, so their lines carry a fixed glass
# reference I0. Brain is not a DAB negative control (it is the highest IL-15 group), so for
# DAB the negative reference is the controls.
#
# DAB_LIST, DAB_OUT and DAB_CHANNEL select the run. The defaults scan DAB over 0.02-0.80 and
# Vector Red over 0.02-0.40, re-measuring the Vector Red scan as a reproducibility check.
# DAB_LIST=dab_endpoint_inputs.txt DAB_OUT=dab_endpoint DAB_CHANNEL=dab scans the other
# analysed sections for scripts/threshold_evidence/dab_threshold_sensitivity.py.
# Mask settings and per-specimen geometry match the cohort, so denominators match.
# Writes only to outputs/threshold_selection/<DAB_OUT>/.
set -uo pipefail
ROOT="${IL15_ROOT:-/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis}"
source "${ROOT}/config/site.env"
MACRO=${ROOT}/macros/il15_threshold_scan.ijm
LIST=${ROOT}/inputs_lists/${DAB_LIST:-dab_check_inputs.txt}
OUTSUB=${DAB_OUT:-dab_check}
CHANNEL=${DAB_CHANNEL:-both}
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST}")
[[ -n "${LINE}" ]] || { echo "ERROR: no input at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
LABEL=${LINE%%|*}
REST=${LINE#*|}
IMG=${REST%%|*}

# optional third field: fixed glass reference for fields with no glass at the corners
I0ARG=""
[[ "${REST}" == *"|"* ]] && I0ARG="|${REST#*|}"
[[ -f "${IMG}" ]] || { echo "ERROR: image not found: ${IMG}" >&2; exit 1; }
OUTDIR=${ROOT}/outputs/threshold_selection/${OUTSUB}/${LABEL}

# Per-specimen sampling geometry, from the same file the production pipeline reads.
# Whole-section values depend on the mask erosion, so small specimens need their own settings.
GRID=1250; ERODE=100; MINF=0.999; FOVS=500
OVR=${ROOT}/data/decisions/fov_sampling_overrides.csv
OLINE=$(awk -F, -v s="${LABEL}" 'NR>1 && $2==s {print $4","$5","$6","$7; exit}' "${OVR}")
if [[ -n "${OLINE}" ]]; then
    IFS=, read -r g e m f <<< "${OLINE}"
    GRID=${g:-$GRID}; ERODE=${e:-$ERODE}; MINF=${m:-$MINF}; FOVS=${f:-$FOVS}
    echo "sampling override: pitch ${GRID} um, erosion ${ERODE} um, min tissue ${MINF}, field ${FOVS} um"
fi
mkdir -p "${OUTDIR}"
echo "task ${SLURM_ARRAY_TASK_ID} | ${LABEL} | node $(hostname) | $(date)"
IJPRIV="${TMPDIR:-/tmp}/ij_${SLURM_ARRAY_TASK_ID:-0}_$$"
mkdir -p "${IJPRIV}/prefs"
export _JAVA_OPTIONS="-Dnet.imagej.legacy.SingleInstance=false -Djava.io.tmpdir=${IJPRIV} -Djava.util.prefs.userRoot=${IJPRIV}/prefs -Djava.util.prefs.systemRoot=${IJPRIV}/prefs"
JAVA_MEM=$(( ${SLURM_MEM_PER_NODE:-49152} / 1024 - 8 ))
rm -f /tmp/ImageJ-${USER}-*.stub 2>/dev/null || true
# -batch, not -macro: a macro error is printed instead of hanging in an invisible dialog.
# Every positional argument carries a real value (ImageJ's split() collapses blanks).
xvfb-run -a "$FIJI" --mem=${JAVA_MEM}g --console -port0 -batch "$MACRO" \
    "${IMG}|${OUTDIR}|${GRID}|${ERODE}|${MINF}|${FOVS}|7|0.10|${CHANNEL}${I0ARG}"
echo "fiji rc=$? $(date)"
ls -la "${OUTDIR}"/*threshold_scan.csv 2>/dev/null | awk '{print $5"  "$NF}'
