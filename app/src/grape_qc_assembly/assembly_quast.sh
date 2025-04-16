#!/bin/bash -u
# for each of strains (SPAdes or Skesa --> quast --> filter contigs by length with seqkit) --> quast summary

QUAST_RESULT_DIR=quast
ASSEMBLY_SUMMARY_FILE=assembly_summary.tsv
CMDNAME=`basename $0`
# ヘルプメッセージを表示する関数
show_help() {
  echo "Usage: $CMDNAME --assembler ASSEMBLER --input-dir INPUT_DIR --assembly-list ASSEMBLY_LIST [--threads THREADS] [--cleanup CLEANUP] [--cutoff-length CUTOFF_LENGTH] [--scaffold SCAFFOLD] [--extension EXT]"
  echo
  echo "Required arguments:"
  echo "  -a, --assembler           Assembler to use (p for SPAdes, k for Skesa)"
  echo "  -i, --input-dir           Input directory"
  echo "  -l, --assembly-list       Path to the assembly fastq list file"
  echo
  echo "Optional arguments:"
  echo "  -c, --cleanup             Cleanup intermediate data like quast raw result [1: cleanup, 0: not cleanup] (Default: 1)"
  echo "  -L, --cutoff-length       Cutoff length for contigs (Default: 500)"
  echo "  -s, --scaffold            Create scaffold by SPAdes [1: create, 0: not create] (Default: 0, if you choose Skesa, this option will be ignored)"
  echo "  -t, --threads             Number of threads (Default: 8)"
  echo "  -x, --extension           Extension for contig files (Default: fasta)"
}

post_process() {
  # quast contig
  local quast_out_dir=$1
  local contig_file=$2

  filtered_contig_file="filtered_${contig_file}"
  source activate seqkit
  seqkit seq -j ${THREADS} -m ${CUTOFF_LENGTH} "${contig_file}" > "${filtered_contig_file}"
  seqkit sort -j ${THREADS} -lr2 "${filtered_contig_file}" > "${contig_file}" # filter + sortしたものでオリジナルを上書き
  conda deactivate
  rm "${filtered_contig_file}" *.fai # remove seqkit information file

  source activate quast
  quast --threads ${THREADS} "${contig_file}" -o "${quast_out_dir}"
  conda deactivate
}

# デフォルト値の設定
THREADS=8
CLEANUP=1
CUTOFF_LENGTH=500
EXT=fasta
SCAFFOLD=0
# 引数の処理
while [[ $# -gt 0 ]]; do
  case "$1" in
    -a|--assembler)
      ASSEMBLER="$2"
      shift 2
      ;;
    -c|--cleanup)
      CLEANUP="$2"
      shift 2
      ;;
    -l|--assembly-list)
      ASSEMBLY_LIST="$2"
      shift 2
      ;;
    -L|--cutoff-length)
      CUTOFF_LENGTH="$2"
      shift 2
      ;;
    -i|--input-dir)
      INPUT_DIR="$2"
      shift 2
      ;;
    -s|--scaffold)
      SCAFFOLD="$2"
      shift 2
      ;;
    -t|--threads)
      THREADS="$2"
      shift 2
      ;;
    -x|--extension)
      EXT="$2"
      shift 2
      ;;
    --help)
      show_help
      ;;
    *)
      echo "Unknown option: $1" >&2
      show_help
      exit 1
      ;;
  esac
done

# 必須引数の確認
if [ -z "$ASSEMBLER" ] || [ -z "$ASSEMBLY_LIST" ] || [ -z "$INPUT_DIR" ]; then
  echo "Error: Missing required arguments!" >&2
  show_help
  exit 1
fi

if [ "$ASSEMBLER" != "P" ] && [ "$ASSEMBLER" != "p" ] && [ "$ASSEMBLER" != "K" ] && [ "$ASSEMBLER" != "k" ]; then
  echo "Error: Invalid assembler option!" >&2
  show_help
  exit 1
fi

# Validate CLEANUP option
if [ $CLEANUP -ne 0 ] && [ $CLEANUP -ne 1 ]; then
  echo "Error: Invalid cleanup option!" >&2
  show_help
  exit 1
fi

# Validate extension option
if [ -z "$EXT" ]; then
  echo "Error: Missing extension option!" >&2
  show_help
  exit 1
fi

# Validate scaffold option
if [ $SCAFFOLD -ne 0 ] && [ $SCAFFOLD -ne 1 ]; then
  echo "Error: Invalid scaffold option!" >&2
  show_help
  exit 1
fi

if [ $SCAFFOLD -eq 1 ] && [ "$ASSEMBLER" != "P" ] && [ "$ASSEMBLER" != "p" ]; then
  echo "Warning: Scaffold option will be ignored because it is valid only for SPAdes."
fi

# オプションの確認（表示などのデバッグ用）
echo "Assembler        : $ASSEMBLER"
echo "Assembly list    : $ASSEMBLY_LIST"
echo "Input directory  : $INPUT_DIR"
echo "Cutoff length    : $CUTOFF_LENGTH"
echo "Cleanup          : $CLEANUP"
echo "Scaffold         : $SCAFFOLD"
echo "Threads          : $THREADS"
echo "Extension        : $EXT"

DIR=${BAC1}/grape_qc_assembly
INPUT_DIR=${INPUT_DIR%/}  # 末尾のスラッシュを削除

mapfile -t lines < "$ASSEMBLY_LIST"
for line in "${lines[@]}"; do
  # 行をタブで分割して配列に格納
  IFS=$'\t' read -r -a columns <<< "$line"
  strain=${columns[0]}
  paired1_file_name=${columns[1]}
  paired2_file_name=${columns[2]}
  unpaired1_file_name=${columns[3]}
  unpaired2_file_name=${columns[4]}

  echo "Assembly: [$strain]"
  # Define input files
  paired1="${INPUT_DIR}/${paired1_file_name}"
  paired2="${INPUT_DIR}/${paired2_file_name}"
  unpaired1=""
  if [ -n "$unpaired1_file_name" ]; then
    unpaired1="${INPUT_DIR}/${unpaired1_file_name}"
  fi
  unpaired2=""
  if [ -n "$unpaired2_file_name" ]; then
    unpaired2="${INPUT_DIR}/${unpaired2_file_name}"
  fi

  # Assembly
  echo "INPUT:"
  echo "  paired1  : $paired1"
  echo "  paired2  : $paired2"
  echo "  unpaired1: $unpaired1"
  echo "  unpaired2: $unpaired2"

  scaffold_file=""
  if [ $ASSEMBLER = "P" ] || [ $ASSEMBLER = "p" ]; then
    echo "SPAdes"
    # SPAdes
    source activate spades
    contig_file="${strain}.${EXT}"
    s_option=("")
    if [ -n "$unpaired1" ]; then
      # SPAdesのunpairedに対する引数は1つなので、unpaired1と2を結合して1つのファイルにする
      merged_unpaired=${strain}_u.fastq.gz
      cat "${unpaired1}" "${unpaired2}" > "${merged_unpaired}"
      s_option=("-s" "${merged_unpaired}")
    fi
    echo "SPAdes: spades.py --careful --cov-cutoff 10 \
      -1 "${paired1}" -2 "${paired2}" -o "${strain}" -t ${THREADS} ${s_option[@]}"
    spades.py --careful --cov-cutoff 10 \
      -1 "${paired1}" -2 "${paired2}" -o "${strain}" -t ${THREADS} ${s_option[@]}
    echo "SPAdes: Done"

    # Extract contig file
    cp ${strain}/contigs.fasta "${contig_file}"
    if [ $SCAFFOLD -eq 1 ]; then
      scaffold_file="${strain}_s.${EXT}"
      cp ${strain}/scaffolds.fasta "${scaffold_file}"
    fi
    rm -r ${strain}
    if [ -n "$merged_unpaired" ]; then
      rm ${merged_unpaired}
    fi
  elif [ $ASSEMBLER = "K" ] || [ $ASSEMBLER = "k" ]; then
    # Skesa
    source activate skesa
    contig_file="${strain}.${EXT}"
    reads_option=("${paired1}" "${paired2}")
    if [ -n "$unpaired1" ]; then
      reads_option+=("${unpaired1}")
    fi
    if [ -n "$unpaired2" ]; then
      reads_option+=("${unpaired2}")
    fi
    skesa --cores ${THREADS} --reads ${reads_option[@]} > "${contig_file}"
  fi
  conda deactivate

  # quast contig
  post_process "${QUAST_RESULT_DIR}/${strain}" "${contig_file}"
  if [ -n "$scaffold_file" ]; then
    post_process "${QUAST_RESULT_DIR}/${strain}_s" "${scaffold_file}"
  fi
done

# quast summary
if [ $ASSEMBLER = "P" ] || [ $ASSEMBLER = "p" ]; then
  # spades summary
  source activate spades
  version_text=$(spades.py --version)
  conda deactivate
elif [ $ASSEMBLER = "K" ] || [ $ASSEMBLER = "k" ]; then
  # skesa summary
  source activate skesa
  version_text=$(skesa --version)
  conda deactivate
fi
source activate original
python "${DIR}/write_assembly_summary.py" "${QUAST_RESULT_DIR}" "${ASSEMBLY_LIST}" "${version_text}" "${ASSEMBLY_SUMMARY_FILE}"
conda deactivate

# remove quast directory
if [ $CLEANUP = "1" ]; then
  rm -r ${QUAST_RESULT_DIR}
fi
