#!/bin/bash
# $1: Input alignment in FASTA format
# $2: Ultrafast bootstrap replicates (>=1000)

set -eo pipefail

CMDNAME=`basename $0`

# parameter check
if [ "$#" -lt 1 ]; then
  echo "Usage: bash $CMDNAME input [bootstrap]"
  echo '  $1  Input alignment'
  echo '  $2  Ultrafast bootstrap replicates (1000 by default)'
  exit 1
fi

if [ ! -f "$1" ]; then
  echo "Error: Input file '$1' does not exist or is not a regular file." >&2
  exit 1
fi

# assign parameter
INPUT_FILE=$1

if [[ -n "$2" ]]; then
  if [[ "$2" =~ ^[1-9][0-9]*$ ]] && [ "$2" -ge 1000 ]; then
    BB=$2
  else
    echo "Error: Bootstrap value parameter '$2' must be an integer >= 1000." >&2
    exit 1
  fi
else
  BB=1000
fi

# remove invariable region
OUT_DIR=iqtree_results_`date +%Y%m%d_%H%M%S`
mkdir "${OUT_DIR}"
INPUT_FILE_NAME=$(basename "${INPUT_FILE}")
VARIABLE_SNP_FILE_NAME="${INPUT_FILE_NAME}"
VARIABLE_SNP_FILE_PATH="${OUT_DIR}/${VARIABLE_SNP_FILE_NAME}"
conda run -n snippy snp-sites -c "${INPUT_FILE}" > "${VARIABLE_SNP_FILE_PATH}"

# iqtree
cd "${OUT_DIR}"
conda run -n iqtree iqtree -s "${VARIABLE_SNP_FILE_NAME}" -nt AUTO -bb "${BB}" -m MFP+ASC

cp -pv "${VARIABLE_SNP_FILE_NAME}".contree "${VARIABLE_SNP_FILE_NAME}".contree.nwk
