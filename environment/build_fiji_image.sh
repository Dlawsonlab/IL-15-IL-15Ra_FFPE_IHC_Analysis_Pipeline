#!/bin/bash
#SBATCH --job-name=il15_fiji_image
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/fiji_image_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis/logs/fiji_image_%j.err

# Record the cluster Fiji install, then build it into an Apptainer image.
#
# Writes, in environment/:
#   fiji_manifest.tsv         SHA-256 and size of every file in Fiji.app
#   fiji_versions.txt         ImageJ jar, Java version, plugin hash check
#   fiji_image.sha256         hash of the built image
#   fiji_image_inspect.txt    the image's labels and definition, as Apptainer reports them
# and the image itself at $FIJI_SIF (default /pub/evaz/containers/).
#
# Submit:  sbatch environment/build_fiji_image.sh

set -euo pipefail

ROOT="${IL15_ROOT:-/share/crsp/lab/dalawson/share/pascal-tim-collab/7_mouse_IL-15Ra_Liver_Brain/IL15_IHC_analysis}"
source "${ROOT}/config/site.env"
ENV="${ROOT}/environment"
FIJI_APP="$(dirname "${FIJI}")"
SIF="${FIJI_SIF:-/pub/${USER}/containers/il15_fiji_$(date +%Y%m%d).sif}"
PLUGIN_SHA=1690e891f5aef7bbe116806f425539c6b00d38e5d9bd7fd11e962026886e1002

echo "Fiji.app: ${FIJI_APP}  ($(du -sh "${FIJI_APP}" | cut -f1))  node $(hostname)  $(date)"

# ---- manifest ----
( cd "${FIJI_APP}"
  printf 'path\tbytes\tsha256\n'
  find . -type f -print0 | sort -z | while IFS= read -r -d '' f; do
      printf '%s\t%s\t%s\n' "${f#./}" "$(stat -c %s "$f")" "$(sha256sum "$f" | cut -d' ' -f1)"
  done ) > "${ENV}/fiji_manifest.tsv"
echo "manifest: $(($(wc -l < "${ENV}/fiji_manifest.tsv") - 1)) files"

# ---- versions ----
{
    echo "# Fiji install recorded $(date +%Y-%m-%d) from ${FIJI_APP}"
    echo "imagej_jar: $(cd "${FIJI_APP}/jars" && ls ij-*.jar)"
    echo "imagej2_jar: $(cd "${FIJI_APP}/jars" && ls imagej-2*.jar 2>/dev/null || echo none)"
    echo "fiji_jar: $(cd "${FIJI_APP}/jars" && ls fiji-2*.jar 2>/dev/null || echo none)"
    echo "bioformats_jar: $(cd "${FIJI_APP}/jars" && ls bio-formats/formats-api-*.jar 2>/dev/null || echo none)"
    JAVA=$(find "${FIJI_APP}/java" -type f -path '*/bin/java' | head -1)
    echo "java: ${JAVA#${FIJI_APP}/}"
    "${JAVA}" -version 2>&1 | sed 's/^/  /'
    CD2=$(find "${FIJI_APP}/plugins" -iname 'colour_deconvolution2*.jar' | head -1)
    echo "colour_deconvolution2: ${CD2#${FIJI_APP}/}"
    echo "  sha256 $(sha256sum "${CD2}" | cut -d' ' -f1)"
    if [[ "$(sha256sum "${CD2}" | cut -d' ' -f1)" == "${PLUGIN_SHA}" ]]; then
        echo "  matches macros/colour_deconvolution2.jar"
    else
        echo "  DIFFERS from macros/colour_deconvolution2.jar"
    fi
} > "${ENV}/fiji_versions.txt"
cat "${ENV}/fiji_versions.txt"

# ---- image ----
module load apptainer/1.4.5
mkdir -p "$(dirname "${SIF}")" "/pub/${USER}/.apptainer_cache"
export APPTAINER_CACHEDIR="/pub/${USER}/.apptainer_cache"
export APPTAINER_TMPDIR="${TMPDIR:-/tmp}"
apptainer --version
apptainer build --fakeroot --force "${SIF}" "${ENV}/fiji.def"

{ echo "# $(basename "${SIF}")  built $(date +%Y-%m-%d) with $(apptainer --version)"
  sha256sum "${SIF}" | sed "s|${SIF}|$(basename "${SIF}")|"; } > "${ENV}/fiji_image.sha256"
apptainer inspect --labels --deffile "${SIF}" > "${ENV}/fiji_image_inspect.txt"
cat "${ENV}/fiji_image.sha256"
ls -lh "${SIF}"
echo "done $(date)"
