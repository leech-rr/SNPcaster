#!/bin/bash

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

FILE_IN=$1
OUT_DIR=raxml_results_`date +%Y%m%d_%H%M%S`
mkdir -v "${OUT_DIR}"

# SIZE
if [[ -n "$3" ]] && [[ "$3" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
  SIZE=$3
else
  SIZE=1000
fi

# run modeltest
source activate raxml-ng
modelttest_prefix="${OUT_DIR}/modeltest-ng"
modeltest-ng -i $FILE_IN -t ml -p $2 --output "${modelttest_prefix}"
modelttest_out="${modelttest_prefix}.out"

DIR_IN=$(dirname ${FILE_IN})
# Extract raxml command with best-fit model from modeltest-ng result(*.out)
cmd=$(sed -n '/raxml-ng/p' "${modelttest_out}" | tail -n1 | sed 's/^[[:space:]]*>[[:space:]]*//')
RAXML_PREFIX="${OUT_DIR}/"
$cmd --prefix "${RAXML_PREFIX}" --all --bs-trees $SIZE --threads auto --workers auto
TREE_FILE="${RAXML_PREFIX}.raxml.support"

conda deactivate

# Change extension to nwk
FILE_FINAL="${OUT_DIR}/$(basename "${FILE_IN}")_bootstrap.nwk"
cp -pv "${TREE_FILE}" "${FILE_FINAL}"
