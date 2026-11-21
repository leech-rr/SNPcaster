#!/bin/bash
# Process 1: modeltest
# Process 2: raxml-ng
# Process 3: change filename
# 
#$1 multi-fasta file
#$2 threads
#$3 optional, bootstrap value (default 1000)

set -eo pipefail

# parameter check
if [ "$#" -lt 2 ]; then
  cat << 'EOF'
Usage: bash raxml-ng.sh $1 $2 $3
  $1 multi-fasta file
  $2 threads
  $3 optional, bootstrap value (default 1000)
EOF
  exit 1
fi

if [ ! -f "$1" ]; then
  echo "Error: Input file '$1' does not exist or is not a regular file." >&2
  exit 1
fi

if ! [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: Threads parameter '$2' must be a positive integer." >&2
  exit 1
fi

INPUT_FILE=$1
# determine SIZE
if [[ -n "$3" ]]; then
  if [[ "$3" =~ ^[1-9][0-9]*$ ]]; then
    SIZE=$3
  else
    echo "Error: Bootstrap value parameter '$3' must be a positive integer (non-zero)." >&2
    exit 1
  fi
else
  SIZE=1000
fi

OUT_DIR=raxml_results_`date +%Y%m%d_%H%M%S`
mkdir -v "${OUT_DIR}"
INPUT_FILE_NAME=$(basename "${INPUT_FILE}")
VARIABLE_SNP_FILE_NAME="${INPUT_FILE_NAME}"
VARIABLE_SNP_FILE_PATH="${OUT_DIR}/${VARIABLE_SNP_FILE_NAME}"
conda run -n snippy snp-sites -c "${INPUT_FILE}" > "${VARIABLE_SNP_FILE_PATH}"


cd "${OUT_DIR}"
source activate raxml-ng
# run modeltest
modeltest_prefix="modeltest-ng"
modeltest-ng -i "${VARIABLE_SNP_FILE_NAME}" -t ml -p "$2" -a lewis -h ug --output "${modeltest_prefix}"
modeltest_out="${modeltest_prefix}.out"

# Extract the RAxML-NG model string selected by BIC
best_model=$(awk '
  /^Best model according to BIC$/ { in_bic = 1; next }
  /^Best model according to /     { in_bic = 0 }
  in_bic && /> raxml-ng / {
    for (i = 1; i <= NF; i++) {
      if ($i == "--model") {
        print $(i + 1)
        exit
      }
    }
  }
' "${modeltest_out}")

if [[ -z "${best_model}" || "${best_model}" != *"+ASC_LEWIS" ]]; then
  echo "Error: Failed to extract a valid BIC model with +ASC_LEWIS from '${modeltest_out}' (got: '${best_model}')." >&2
  exit 1
fi

RAXML_PREFIX="${VARIABLE_SNP_FILE_NAME}"
raxml-ng \
  --msa "${VARIABLE_SNP_FILE_NAME}" \
  --model "${best_model}" \
  --prefix "${RAXML_PREFIX}" \
  --all \
  --bs-trees "${SIZE}" \
  --threads "auto{$2}" \
  --workers auto

# Change extension to nwk
TREE_FILE="${RAXML_PREFIX}.raxml.support"
FILE_FINAL="${RAXML_PREFIX}_bootstrap.nwk"
cp -pv "${TREE_FILE}" "${FILE_FINAL}"
