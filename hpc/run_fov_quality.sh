#!/bin/bash
#SBATCH --job-name=il15_fovqual
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/fovqual_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/fovqual_%j.err
#
# Native-resolution FOV quality scoring (tear / fold / dried).
#
# The job runs on the cluster because the scorer reads every native crop (~1.5 GB).
# Reading that much across the CRSP SMB mount can stall the mount badly enough that even
# `posix_spawn /bin/bash` fails with EIO. Compute nodes read CRSP directly and avoid that
# failure.
#
# The scorer needs only numpy, scipy and PIL. tifffile is absent on the compute nodes,
# which is why the crops are PNGs.
set -uo pipefail
ROOT="${IL15_ROOT:-/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis}"
source "${ROOT}/config/site.env"
cd "${ROOT}" || exit 1

module load python/3.10.2 2>/dev/null || true
python3 -c 'import numpy, scipy; from PIL import Image' || { echo "missing python deps" >&2; exit 1; }

echo "host $(hostname) | $(date)"
echo "crops on disk: $(find per_sample_out -path '*/fov/raw/*_raw.png' | wc -l)"

JOB_START_EPOCH=$(date +%s)
python3 scripts/pipeline/fov_quality_score.py
RC=$?
echo "scorer rc=${RC} $(date)"

# Success requires fresh output. A failed re-run leaves the previous CSV in place, so a
# presence check alone would report OK.
OUT="${ROOT}/outputs/tables/fov_quality_flagged.csv"
if [[ -f "${OUT}" && -n "$(find "${OUT}" -newermt "@${JOB_START_EPOCH}" 2>/dev/null)" ]]; then
    echo "OK: ${OUT} written fresh"
    wc -l "${ROOT}/outputs/tables/fov_quality_scores.csv" "${OUT}"
else
    echo "FAILED: no fresh ${OUT} (rc=${RC})" >&2; exit 1
fi
