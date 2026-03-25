#!/bin/bash
# fastpによるトリミングを実行するスクリプト (純粋なfastp呼び出しのみ)
# Arguments:
# $1: Input Directory
# $2: Strain List Path
# $3: Fastp Output Directory (fastq files)
# $4: Threads for fastp (3 by default)

CMDNAME=`basename $0`

# 引数チェック
if [ "$#" -lt 3 ]; then
  echo "Usage: bash $CMDNAME input_dir/ list output/fastp [threads]"
  echo '  $1  Input Directory'
  echo '  $2  Strain List Path'
  echo '  $3  Fastp Output Directory (fastq files)'
  echo '  $4  Threads for fastp (3 by default)'
  exit 1
fi

# variables
DIR=${BAC1}/grape_qc_assembly

# 引数代入
INPUT_DIR=$1
LIST=$2
FASTP_DIR=$3
if [ -n "$FASTP_DIR" ]; then
  	FASTP_DIR="${FASTP_DIR%/}"
fi

THREADS=3
if [ -n "$4" ]; then
  THREADS=$4
fi

echo `date +%Y%m%d_%H%M%S`

mapfile -t lines < "$LIST"
for line in "${lines[@]}"; do
  # 行をタブで分割して配列に格納
  IFS=$'\t' read -r -a columns <<< "$line"
  strain=${columns[0]}
  read1_file_name=${columns[1]}
  read2_file_name=${columns[2]}

  in1=${INPUT_DIR}/${read1_file_name}
  in2=${INPUT_DIR}/${read2_file_name}
  out1=${FASTP_DIR}/${strain}_1.fastq.gz
  out2=${FASTP_DIR}/${strain}_2.fastq.gz
  unpaired1=${FASTP_DIR}/${strain}_u1.fastq.gz
  unpaired2=${FASTP_DIR}/${strain}_u2.fastq.gz
  
  source activate fastp
  fastp --thread ${THREADS} --in1 "${in1}" --in2 "${in2}" \
      --out1 "${out1}" --out2 "${out2}" \
      --unpaired1 "${unpaired1}" --unpaired2 "${unpaired2}" \
      -q 30
  conda deactivate
  
  rm -f fastp.json fastp.html 
done

echo `date +%Y%m%d_%H%M%S`
