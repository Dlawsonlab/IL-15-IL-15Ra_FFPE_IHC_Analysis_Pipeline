#!/bin/bash
#SBATCH --job-name=il15_thrscan
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:50:00
#SBATCH --output=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/thrscan_%A_%a.out
#SBATCH --error=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/thrscan_%A_%a.err
#
# Magenta threshold scan against the measured stain vectors.
#
# The measured single-stain vectors (D.M = 0.884) remove most of the DAB-into-magenta
# bleed. A threshold tuned to suppress that bleed would discard faint genuine signal, so
# the threshold is chosen from this scan.
#
# The scan covers 0.02-0.40 in 0.02 steps. The brain floor is the decisive evidence. Brain
# carries no IL-15Ra, so whatever it reads at a given threshold is residual crosstalk. The
# defensible choice is the lowest threshold that keeps brain at the floor.
#
# The job writes to a separate tree (thrscan/), so no cohort output is touched.
set -uo pipefail
ROOT="${IL15_ROOT:-/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis}"
source "${ROOT}/config/site.env"
MACRO=${ROOT}/macros/il15_threshold_scan.ijm

IDX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ROOT}/inputs_lists/thrscan_samples.txt")
[[ -n "${IDX}" ]] || { echo "ERROR: no sample at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IMG=$(sed -n "${IDX}p" "${ROOT}/imagelist.txt")
SUB=$(sed -n "${IDX}p" "${ROOT}/subdirs.txt")
[[ -f "${IMG}" ]] || { echo "ERROR: no image at line ${IDX}" >&2; exit 1; }

OUTDIR=${ROOT}/thrscan/${SUB}
mkdir -p "${OUTDIR}"
echo "task ${SLURM_ARRAY_TASK_ID} | ${SUB} | node $(hostname) | $(date)"

IJPRIV="${TMPDIR:-/tmp}/ij_${SLURM_ARRAY_TASK_ID:-0}_$$"
mkdir -p "${IJPRIV}/prefs"
export _JAVA_OPTIONS="-Dnet.imagej.legacy.SingleInstance=false -Djava.io.tmpdir=${IJPRIV} -Djava.util.prefs.userRoot=${IJPRIV}/prefs -Djava.util.prefs.systemRoot=${IJPRIV}/prefs"
JAVA_MEM=$(( ${SLURM_MEM_PER_NODE:-49152} / 1024 - 8 ))

rm -f /tmp/ImageJ-${USER}-*.stub 2>/dev/null || true
# The mask settings match the cohort (chroma 7, OD 0.10), so the denominator matches the
# cohort and only the threshold varies. Every positional argument carries a real value,
# because ImageJ's split() collapses consecutive delimiters and would shift them silently.
xvfb-run -a "$FIJI" --mem=${JAVA_MEM}g --console -port0 -batch "$MACRO" \
    "${IMG}|${OUTDIR}|1250|100|0.999|500|7|0.10"
echo "fiji rc=$? $(date)"
ls -la "${OUTDIR}"/*_threshold_scan.csv 2>/dev/null | awk '{print $5"  "$NF}'
