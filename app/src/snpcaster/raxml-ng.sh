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
modeltest-ng -i "${VARIABLE_SNP_FILE_NAME}" -t ml -p $2 --output "${modeltest_prefix}"
modeltest_out="${modeltest_prefix}.out"
# Extract raxml command from ModelTest-NG output
raw_cmd=$(sed -n '/raxml-ng/p' "${modeltest_out}" | tail -n1 | sed 's/.*>[[:space:]]*//')
# add +ASC_LEWIS to --model option to optimize SNP phylogenetic tree
cmd=$(echo "$raw_cmd" | sed 's/--model \(\S*\)/--model \1+ASC_LEWIS/')
RAXML_PREFIX="${VARIABLE_SNP_FILE_NAME}"
$cmd --prefix "${RAXML_PREFIX}" --all --bs-trees $SIZE --threads auto --workers auto

# Change extension to nwk
TREE_FILE="${RAXML_PREFIX}.raxml.support"
FILE_FINAL="${RAXML_PREFIX}_bootstrap.nwk"
cp -pv "${TREE_FILE}" "${FILE_FINAL}"
