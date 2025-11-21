#!/bin/bash
# bactsnp, snps.vcf, snippy, gubbins
set -euo pipefail

CMDNAME=$(basename $0)
echo "${CMDNAME}"
DIR0="${BAC1}"/snpcaster
source "${DIR0}"/snpcaster_base.sh

function _usage() {
  cat << __EOF__

  Usage: $(basename "$0") [OPTION]...

  -i      Strain List (Text file listing strain names)   (Required!)
  -r      Reference Sequence File                        (Required!)
  -a      Allele frequency threshold for BactSNP. If allele frequency is smaller than this value, an ambiguous allele is called. (Default: 0.9)
  -c      Gap threshold for clustered SNPs. If set to 0, no SNPs will be removed. (Default: 0)
  -d      Mask file name
  -f      Fastq List (Tab-separated file with strain name, read1, and read2 file names) (Optional, but required for BactSNP execution)
  -g      Run Gubbins [1: run, 0: don't run]              (Default: 0)
  -j      Number of Jobs for BactSNP                      (Default: 4)
  -t      Number of Threads                               (Default: 8)

  [File Format Specifications]

  * Strain List (-i option): A text file listing strain names, one per line.
    Example:
      strain1
      strain2
      ...

  * Fastq List (-f option): A tab-separated file with strain name, read1 file name, and read2 file name.
    Example:
      strain1  strain1_R1.fastq.gz  strain1_R2.fastq.gz
      strain2  strain2_R1.fastq.gz  strain2_R2.fastq.gz
      ...

  (Note: Ensure that file extensions for Fastq files are correct.)

__EOF__
}

function create_report_header() {
  local REPORT=$1
  echo "This is SNPcaster " > "${REPORT}"
  echo "Version : $(get_version)" >> "${REPORT}"
  echo "Strain list : ${LIST_FULL}" >> "${REPORT}"
  echo "Fastq list : ${FASTQ_LIST}" >> "${REPORT}"
  echo "Mask file : ${MASK_FILE}" >> "${REPORT}"
  echo "Reference file : ${REF}" >> "${REPORT}"
  echo "Clustered SNP removal (bp) : ${GAP}" >> "${REPORT}"
  echo "Threads : ${THREAD}" >> "${REPORT}"
  echo "BactSNP allele-freq : ${ALLELE_FREQ}" >> "${REPORT}"
  echo "BactSNP jobs : ${JOBS}" >> "${REPORT}"
  if [ $EXEC_GUBBINS -eq 1 ]; then
    echo "Recombinogenic region detection (Gubbins) : ON" >> "${REPORT}"
  else
    echo "Recombinogenic region detection (Gubbins) : OFF" >> "${REPORT}"
  fi
}

DIR0=$BAC1/snpcaster

REF=""
ALLELE_FREQ=0.9
LIST_FULL=""
GAP=0
MASK_FILE=""
THREAD=8
JOBS=4
EXEC_GUBBINS=0
FASTQ_LIST=""
while getopts "d:i:r:a:c:f:g:j:t:" optKey; do
  case "$optKey" in
  d)
    MASK_FILE=${OPTARG}
    ;;
  i)
    LIST_FULL=${OPTARG}
    ;;
  r)
    REF=${OPTARG}
    ;;
  a)
    ALLELE_FREQ=${OPTARG}
    ;;
  c)
    GAP=${OPTARG}
    ;;
  f)
    FASTQ_LIST=${OPTARG}
    ;;
  g)
    EXEC_GUBBINS=${OPTARG}
    ;;
  j)
    JOBS=${OPTARG}
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

if [ ! -r "${REF}" ]; then
  echo "Can not find or read reference file [${REF}]." >&2
  _usage
  exit 1
fi

if [ ! -r "${LIST_FULL}" ]; then
  echo "Can not find or read strain list file [${LIST_FULL}]" >&2
  _usage
  exit 1
fi
using_mask_file=0
if [ -n "${MASK_FILE:-}" ]; then
  if [ ! -r "${MASK_FILE}" ]; then
    echo "Can not find or read mask file [${MASK_FILE}]" >&2
    _usage
    exit 1
  fi
  using_mask_file=1
fi

if [ -n "${FASTQ_LIST:-}" ] && [ ! -r "${FASTQ_LIST}" ]; then
  echo "Error: Can not find or read fastq list [${FASTQ_LIST}]." >&2
  _usage
  exit 1
fi

LIST_FILE_NAME=${LIST_FULL##*/}
OUTPUT_DIR=snpcaster_$(date +%Y%m%d_%H%M%S)_${LIST_FILE_NAME}
echo "MASK FILE: $MASK_FILE"
echo "STRAIN LIST: $LIST_FULL"
echo "REFERENCE FILE: $REF"
echo "ALLELE FREQUENCY: $ALLELE_FREQ"
echo "CLUSTER GAP: $GAP"
echo "FASTQ LIST: $FASTQ_LIST"
echo "THREAD: $THREAD"
echo "GUBBINS FLG: $EXEC_GUBBINS"
echo "JOBS: $JOBS"
echo "OUTPUT_DIR: $OUTPUT_DIR"

# Create an output folder
[ -d "${OUTPUT_DIR}" ] && rm -r "${OUTPUT_DIR}"
mkdir -v "${OUTPUT_DIR}"

BACTSNP_LIST=bactsnp_list
BACTSNP_SKIP_LIST=bactsnp_skip_list
prepare4snpcaster "${BACTSNP_LIST}" "${FASTQ_LIST}" "${LIST_FULL}" "${OUTPUT_DIR}" "${BACTSNP_SKIP_LIST}" "${REF}" "${MASK_FILE}"

cd "${OUTPUT_DIR}" || exit 1
pwd

REF_FILE_NAME=$(basename "${REF}")
LIST_FILE_NAME=$(basename "${LIST_FULL}")
MASK_FILE_NAME=$(basename "${MASK_FILE}")

process_bactsnp "${BACTSNP_LIST}" "${REF_FILE_NAME}" "${ALLELE_FREQ}" "${THREAD}" "${JOBS}"
process_snippy "${LIST_FILE_NAME}" "${REF_FILE_NAME}"
postprocess_snippy "${GAP}" "${THREAD}" "${MASK_FILE_NAME}"
if [ "${EXEC_GUBBINS}" -eq 1 ]; then
  process_gubbins "${THREAD}"
  postprocess_gubbins "${THREAD}"
fi

REPORT=SNPcaster_report.txt
create_report_header "${REPORT}"
append_report_information "${REPORT}" 1 "${GAP}" "${using_mask_file}" "${EXEC_GUBBINS}"
# remove unnecessary data
if [ -f "${BACTSNP_SKIP_LIST}" ]; then
  for strain in $(cat "${BACTSNP_SKIP_LIST}"); do
    rm -r "${strain:?}"/
  done
fi
rm -v matching_list.txt core.ref.fa

# Rename fasta file extension to ${EXT}
echo '=================== Rename .fa and .fas to .fasta ==================='
if [ -d "${BACTSNP_OUTDIR}" ]; then
  rename_fasta_extension_in_directories "${BACTSNP_OUTDIR}"
fi

rename_fasta_extension_in_directories "${SNIPPY_OUTDIR}" "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"

if [ "${EXEC_GUBBINS}" -eq 1 ]; then
  rename_fasta_extension_in_directories "${GUBBINS_OUTDIR}" "${RESULTS_WITH_GUBBINS_OUTDIR}"
fi

echo "Done!"
