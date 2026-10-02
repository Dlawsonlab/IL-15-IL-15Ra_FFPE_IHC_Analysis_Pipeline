#!/bin/bash
#SBATCH --job-name=il15_fiji_image_test
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/fiji_image_test_%A_%a.out
#SBATCH --error=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/fiji_image_test_%A_%a.err
#SBATCH --array=21,59

# Rerun the whole-slide macro on one section inside the Fiji image, then compare its CSVs
# with the cohort run in per_sample_out/. Identical CSVs show that the image reproduces the
# reported numbers. The array index is the line of imagelist.txt / subdirs.txt, as in
# hpc/run_il15_array.sh; the defaults are two small sections (a BCBM, a normal liver).
#
# Submit:  sbatch environment/test_fiji_image.sh
#          FIJI_SIF=/path/to/image.sif sbatch --array=<N> environment/test_fiji_image.sh
#
# Writes only under /pub/$USER/il15_fiji_image_test/. Never touches per_sample_out/.

set -uo pipefail

ROOT="${IL15_ROOT:-/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis}"
SIF="${FIJI_SIF:-$(ls -t /pub/${USER}/containers/il15_fiji_*.sif | head -1)}"
MACRO="${ROOT}/macros/il15_ihc_pipeline.ijm"
N="${SLURM_ARRAY_TASK_ID}"
IMG=$(sed -n "${N}p" "${ROOT}/imagelist.txt")
SUB=$(sed -n "${N}p" "${ROOT}/subdirs.txt")
REF="${ROOT}/per_sample_out/${SUB}"
OUT="/pub/${USER}/il15_fiji_image_test/${SUB}/"
rm -rf "${OUT}"; mkdir -p "${OUT}"

# Sampling overrides, read as hpc/run_il15_array.sh reads them.
EXTRA=""
LINE=$(awk -F, -v s="${SUB}" 'NR>1 && $2==s {print $4","$5","$6","$7; exit}' \
    "${ROOT}/data/decisions/fov_sampling_overrides.csv")
if [[ -n "${LINE}" ]]; then
    EXTRA="|$(echo "${LINE}" | tr ',' '|')"
fi

echo "task ${N} | ${SUB} | image $(basename "${SIF}") | node $(hostname) | $(date)"

module load apptainer/1.4.5
# Same launch as the cohort: HOME at the project root, a private Java prefs root, no
# single-instance handoff, -port0, xvfb and never --headless.
IJPRIV="${TMPDIR:-/tmp}/ij_${N}_$$"
mkdir -p "${IJPRIV}/prefs"
JAVA_MEM="$(( ${SLURM_MEM_PER_NODE:-32000} / 1000 - 6 ))g"
apptainer exec --cleanenv --bind /share/crsp,/pub/${USER},"${IJPRIV}" \
    --home "${ROOT}" \
    --env _JAVA_OPTIONS="-Dnet.imagej.legacy.SingleInstance=false -Djava.io.tmpdir=${IJPRIV} -Djava.util.prefs.userRoot=${IJPRIV}/prefs -Djava.util.prefs.systemRoot=${IJPRIV}/prefs" \
    "${SIF}" xvfb-run -a /opt/Fiji.app/ImageJ-linux64 --mem="${JAVA_MEM}" --console -port0 \
    -macro "${MACRO}" "${IMG}|${OUT}${EXTRA}"
echo "launcher rc=$?"

# ---- compare with the cohort run ----
status=0
for ref in "${REF}"/*_absorbance.csv "${REF}"/fov/*_fov.csv; do
    new="${OUT}${ref#${REF}/}"
    if [[ ! -f "${new}" ]]; then
        echo "MISSING  ${ref#${REF}/}"; status=1
    elif cmp -s "${ref}" "${new}"; then
        echo "IDENTICAL  ${ref#${REF}/}"
    else
        echo "DIFFERS  ${ref#${REF}/}"; diff "${ref}" "${new}" | head -20; status=1
    fi
done
echo "done $(date)"
exit ${status}
