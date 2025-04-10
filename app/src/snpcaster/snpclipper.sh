#!/bin/bash
# bactsnp, snps.vcf, snippy, gubbins

CMDNAME=$(basename $0)
echo "${CMDNAME}"
DIR0="${BAC1}"/snpcaster
source "${DIR0}"/snpcaster_base.sh

function _usage() {
  cat << __EOF__

  Usage: $(basename "$0") [OPTION]...

  -i      SNPcaster directory having SNPcaster outputs   (Required!)
  -c      Gap threshold for cluster SNPs. If set to 0, no SNPs will be removed. (Default: 0)
  -d      Mask file name
  -t      Number of Threads                              (Default: 8)

__EOF__
}

function create_report_header() {
  local REPORT=$1
  echo "This is SNPclipper " > "${REPORT}"
  echo "SNPcaster directory : ${SNPCASTER_DIR}" >> "${REPORT}"
  echo "Mask file : ${MASK_FILE}" >> "${REPORT}"
  echo "Clustered SNP removal (bp) : ${GAP}" >> "${REPORT}"
  echo "Threads : ${THREAD}" >> "${REPORT}"
  if [ $EXEC_GUBBINS -eq 1 ]; then
    echo "Recombinogenic region detection (Gubbins) : ON" >> "${REPORT}"
  else
    echo "Recombinogenic region detection (Gubbins) : OFF" >> "${REPORT}"
  fi
}

DIR0=$BAC1/snpcaster

MASK_FILE=""
SNPCASTER_DIR=""
THREAD=8
GAP=0
while getopts "d:i:c:t:" optKey; do
  case "$optKey" in
  d)
    MASK_FILE=${OPTARG}
    ;;
  i)
    SNPCASTER_DIR=${OPTARG}
    ;;
  c)
    GAP=${OPTARG}
    ;;
  t)
    THREAD=${OPTARG}
    ;;
  *) #該当なし
    _usage
    exit 1
    ;;
  esac
done

if [ -z "${SNPCASTER_DIR}" ]; then
  echo "-i option for SNPcaster directory is required." >&2
  _usage
  exit 1
fi

if [ ! -d "${SNPCASTER_DIR}" ]; then
  echo "Can not find SNPcaster output dir [${SNPCASTER_DIR}]" >&2
  _usage
  exit 1
fi
using_mask_file=0
if [ -n "${MASK_FILE}" ]; then
  if [ ! -r "${MASK_FILE}" ]; then
    echo "Can not find or read mask file [${MASK_FILE}]" >&2
    _usage
    exit 1
  fi
  using_mask_file=1
fi

SNPCASTER_DIR_FULL_PATH=$(realpath "${SNPCASTER_DIR}")
OUTPUT_DIR="${SNPCASTER_DIR_FULL_PATH}"/snpclipper_$(date +%Y%m%d_%H%M%S)
# overwrite SNPcaster output directory defined in snpcaster_base.sh
# These folders are input for this script
BACTSNP_OUTDIR="${SNPCASTER_DIR_FULL_PATH}"/"${BACTSNP_OUTDIR}"
SNIPPY_OUTDIR="${SNPCASTER_DIR_FULL_PATH}"/"${SNIPPY_OUTDIR}"
GUBBINS_OUTDIR="${SNPCASTER_DIR_FULL_PATH}"/"${GUBBINS_OUTDIR}"
EXEC_GUBBINS=0
if [ -d "${GUBBINS_OUTDIR}" ]; then
  EXEC_GUBBINS=1
fi

echo "SNPCASTER FOLDER: $SNPCASTER_DIR"
echo "MASK FILE: $MASK_FILE"
echo "CLUSTER GAP: $GAP"
echo "THREAD: $THREAD"
echo "OUTPUT_DIR: $OUTPUT_DIR"

# Create an output folder
[ -d "${OUTPUT_DIR}" ] && rm -r "${OUTPUT_DIR}"
mkdir -v "${OUTPUT_DIR}"
if [ "${using_mask_file}" -eq 1 ]; then
  cp -pv "${MASK_FILE}" "${OUTPUT_DIR}"
fi

cd "${OUTPUT_DIR}" || exit 1
MASK_FILE_NAME=$(basename "${MASK_FILE}")

postprocess_snippy "${GAP}" "${THREAD}" "${MASK_FILE_NAME}"
if [ "${EXEC_GUBBINS}" -eq 1 ]; then
  postprocess_gubbins "${THREAD}"
fi

REPORT=snpclipper_report.txt
create_report_header "${REPORT}"
append_report_information "${REPORT}" 1 "${GAP}" "${using_mask_file}" "${EXEC_GUBBINS}"

rm -v matching_list.txt

# Rename fasta file extension to ${EXT}
echo '===================Rename .fa and .fas to .fasta ==================='
rename_fasta_extension_in_directories "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"
if [ "${EXEC_GUBBINS}" -eq 1 ]; then
  rename_fasta_extension_in_directories "${RESULTS_WITH_GUBBINS_OUTDIR}"
fi

echo "Done!"
