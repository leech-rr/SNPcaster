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
  -t      Number of Threads        (Default: 8)

__EOF__
}

function create_report_header() {
  local REPORT=$1
  echo "This is Gubbins result in SNPcaster " > "${REPORT}"
  source activate gubbins
  echo "Gubbins version : $(run_gubbins.py --version)" >> "${REPORT}"
  conda deactivate
  echo "SNPcaster directory : ${SNPCASTER_DIR}" >> "${REPORT}"
  echo "Threads : ${THREAD}" >> "${REPORT}"
}

DIR0=$BAC1/snpcaster

SNPCASTER_DIR=""
THREAD=8
while getopts "i:t:" optKey; do
  case "$optKey" in
  i)
    SNPCASTER_DIR=${OPTARG}
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

echo "SNPCASTER FOLDER: $SNPCASTER_DIR"
echo "THREAD: $THREAD"
cd "${SNPCASTER_DIR}" || exit 1

# Check if Gubbins results already exist
if [ -d "${GUBBINS_OUTDIR}" ] || [ -d "${RESULTS_WITH_GUBBINS_OUTDIR}" ]; then
  if [ -d "${GUBBINS_OUTDIR}" ]; then
    echo "Gubbins folder [${GUBBINS_OUTDIR}] already exists." >&2
  elif [ -d "${RESULTS_WITH_GUBBINS_OUTDIR}" ]; then
    echo "Results with Gubbins folder [${RESULTS_WITH_GUBBINS_OUTDIR}] already exists." >&2
  fi
  echo "Please remove it, then run this script again to re-run Gubbins." >&2
  echo "Or, use snpclipper.sh to process the Gubbins results by applying filters for the cluster SNP and/or masking regions." >&2
  _usage
  exit 1
fi

process_gubbins "${THREAD}"
postprocess_gubbins "${THREAD}"

REPORT=gubbins_report.txt
create_report_header "${REPORT}"
append_report_information "${REPORT}" 1 0 0 1

# Rename fasta file extension to ${EXT}
echo '===================Rename .fa and .fas to .fasta ==================='
rename_fasta_extension_in_directories "${GUBBINS_OUTDIR}" "${RESULTS_WITH_GUBBINS_OUTDIR}"

echo "Done!"
