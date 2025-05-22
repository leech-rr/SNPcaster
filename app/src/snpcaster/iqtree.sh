#!/bin/bash
# $1: Input alignment in PHYLIP/FASTA/NEXUS/CLUSTAL/MSF format
# $2: Ultrafast bootstrap replicates (>=1000)

CMDNAME=`basename $0`

# parameter check
if [ "$#" -lt 1 ]; then
  echo "Usage: bash $CMDNAME input (100)"
  echo '  $1  Input alignment'
  echo '  $2  Ultrafast bootstrap replicates (1000 by default)'
  exit 1
fi

# assign parameter
INPUT=$1
if [ -n "$2" ] && [ "$2" -ne 0 ]; then
  BB=$2
else
  BB=1000
fi

# iqtree
source activate iqtree
iqtree -s $INPUT -nt AUTO -bb $BB
conda deactivate

OUT=iqtree_results_`date +%Y%m%d_%H%M%S`
mkdir "${OUT}"
mv "${INPUT}".* "${OUT}"/

INPUT_FILENAME=$(basename "${INPUT}")
cp -pv "${OUT}"/"${INPUT_FILENAME}".contree "${OUT}"/"${INPUT_FILENAME}".contree.nwk
