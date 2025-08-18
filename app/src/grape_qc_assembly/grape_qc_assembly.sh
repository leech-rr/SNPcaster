#!/bin/bash -eu

CMDNAME=`basename $0`
echo $CMDNAME

DIR=${BAC1}/grape_qc_assembly
source "${DIR}/strain_list_functions.sh"

function _usage(){
	cat <<__EOF__

  Usage: $(basename "$0") [OPTION]...
  -i Strain short read List 				(Required!)
  -a assembler [p:spades, k:skesa] ( p by default )
  -c configuration file 		
  -L Cutoff length for contigs (Default: 500)
  -p fastp [1:run, 0:not run] ( 1 by default )
  -r reduce memory usage (execute checkm with --reduced-tree option) [1:reduce, 0:not reduce] ( 0 by default )
  -s Create scaffold by SPAdes [1: create, 0: not create] (0 by default, if you choose Skesa, this option will be ignored)
  -t No. of Thread ( 8 by default )
  -u Cleanup intermediate data like quast raw result [1: cleanup, 0: not cleanup] ( 1 by default )
  
__EOF__
}	

# Assembly+Quast and Checkm
# This function needs several variables to be set before calling
# - FASTP: Run fastp or not (1: run, 0: not run)
# - FASTP_DIR: Fastp result directory
# - INITIAL_DIR: Initial directory
# - ASSEMBLER: Assembler to use (p for SPAdes, k for Skesa)
# - THREAD: Number of threads
# - OVER_COVERAGE_STRAIN_LIST: list of strains over coverage
# - EXT: Extension for contig files (Default: fasta)
# - SCAFFOLD: Create scaffold by SPAdes [1: create, 0: not create] (0 by default, if you choose Skesa, this option will be ignored)
# - CLEANUP: Cleanup intermediate data like quast raw result [1: cleanup, 0: not cleanup] (1 by default)
# - CUTOFF_LENGTH: Cutoff length for contigs (Default: 500)
# - REDUCE_MEMORY: Reduce memory usage (execute checkm with --reduced-tree option) [1: reduce, 0: not reduce]
# - OUTPUT_DIR: Output directory
# Following variables will be set after calling this function
# - CHECKM_RESULT: CheckM result file
function assembly_and_qc() {
  # Assembly + Quast
  ASSEMBLY_LIST="list_assembly.tsv"
  > "${ASSEMBLY_LIST}"
  if [ $FASTP = "1" ]; then
    input_dir=${FASTP_DIR}
    for strain in `cat "$OVER_COVERAGE_STRAIN_LIST"`; do
      echo -e "${strain}\t${strain}_1.fastq.gz\t${strain}_2.fastq.gz\t${strain}_u1.fastq.gz\t${strain}_u2.fastq.gz" >> "${ASSEMBLY_LIST}"
    done
  else
    input_dir=${INITIAL_DIR}
    filter_fastq_list "${LIST}" "${OVER_COVERAGE_STRAIN_LIST}" > "${ASSEMBLY_LIST}"
  fi
  bash ${DIR}/assembly_quast.sh -a ${ASSEMBLER} -t ${THREAD} -i "${input_dir}" -l "${ASSEMBLY_LIST}" -x ${EXT} -s ${SCAFFOLD} -c ${CLEANUP} -L ${CUTOFF_LENGTH}

  # checkM
  CHECKM_DIR="${OUTPUT_DIR}/checkm"
  mkdir "${CHECKM_DIR}"
  CHECKM_INPUT_DIR="${CHECKM_DIR}/input"
  CHECKM_OUTPUT_DIR="${CHECKM_DIR}"
  CHECKM_RESULT="${CHECKM_OUTPUT_DIR}/checkm_results.txt"
  # preparation
  if [ ! -d ${CHECKM_INPUT_DIR} ]; then
    mkdir ${CHECKM_INPUT_DIR}
  else
    echo "Error: ${CHECKM_INPUT_DIR} folder already exists!"
    exit 1
  fi

  # copy contig files to checkm folder
  for strain in `extract_strain_names "${ASSEMBLY_LIST}"`; do
    ln -vs "$(pwd)/${strain}.${EXT}" "${CHECKM_INPUT_DIR}/${strain}.${EXT}"
  done
  reduced_tree_option=()
  if [ $REDUCE_MEMORY = "1" ]; then
    reduced_tree_option=("--reduced_tree")
  fi
  source activate checkm
  checkm lineage_wf -t $THREAD -x ${EXT} "${CHECKM_INPUT_DIR}" "${CHECKM_OUTPUT_DIR}" -f ${CHECKM_RESULT} --tab_table ${reduced_tree_option[@]}
  conda deactivate
}

LIST=""
CONFIG=""
THREAD="8"
FASTP="1"
ASSEMBLER="p"
REDUCE_MEMORY="0"
SCAFFOLD="0"
CUTOFF_LENGTH="500"
CLEANUP="1"
while getopts "a:c:i:L:p:r:t:u:s:" optKey; do
  case "$optKey" in
    a)
     ASSEMBLER=${OPTARG}
     ;;
    c)
     CONFIG=${OPTARG}
     ;;
    i)
     LIST=${OPTARG}
     ;;
    L)
     CUTOFF_LENGTH=${OPTARG}
     ;;
    p)
     FASTP=${OPTARG}
     ;;
    r)
     REDUCE_MEMORY=${OPTARG}
     ;;
    s)
     SCAFFOLD=${OPTARG}
     ;;
    t)
     THREAD=${OPTARG}
     ;; 
    u)
     CLEANUP=${OPTARG}
     ;;
    *)  #該当なし
     _usage
     exit 1
     ;;
  esac
done


if [ -z $LIST ]; then
  	_usage
  	exit 0
fi

if [ "$ASSEMBLER" != "P" ] && [ "$ASSEMBLER" != "p" ] && [ "$ASSEMBLER" != "K" ] && \
[ "$ASSEMBLER" != "k" ]; then
	echo 'Enter p or k for ASSEMBLER'
	exit 1
fi

if [ "$FASTP" != "1" ] && [ "$FASTP" != "0" ]; then
	echo 'Enter 1 or 0 for FASTP'
	exit 1
fi

if [ "$REDUCE_MEMORY" != "1" ] && [ "$REDUCE_MEMORY" != "0" ]; then
	echo 'Enter 1 or 0 for REDUCE_MEMORY'
	exit 1
fi

if [ $(($THREAD*1)) -lt 1 ]; then   # 数値かどうかの判定
	echo 'Enter Integer for THREAD'
	exit 1
fi

# validate CLEANUP option
if [ $CLEANUP -ne 0 ] && [ $CLEANUP -ne 1 ]; then
  echo "Error: Invalid cleanup option." >&2
  exit 1
fi

# validate SCAFFOLD option
if [ $SCAFFOLD -ne 0 ] && [ $SCAFFOLD -ne 1 ]; then
  echo "Error: Invalid scaffold option." >&2
  exit 1
fi


# variables
OVER_COVERAGE_STRAIN_LIST=list_over_coverage
rm -f $OVER_COVERAGE_STRAIN_LIST

INITIAL_DIR=$(pwd)
EXT=fasta
ASSEMBLY_SUMMARY_FILE=assembly_summary.tsv

OUTPUT_DIR="${INITIAL_DIR}/qc_assembly_`date +%Y%m%d_%H%M%S`_${LIST}"

echo "LIST:" $LIST
echo "CONFIG:" $CONFIG
echo "THREAD:" $THREAD
echo "FASTP:" $FASTP
echo "ASSEMBLER:" $ASSEMBLER
echo "REDUCE_MEMORY:" $REDUCE_MEMORY
echo "CUTOFF_LENGTH:" $CUTOFF_LENGTH
echo "SCAFFOLD:" $SCAFFOLD
echo "CLEANUP:" $CLEANUP
echo "OUTPUT_DIR:" $OUTPUT_DIR

# Parse configuration file
qc_summary_options=()
if [ -n "$CONFIG" ]; then
  echo "CONFIG=$CONFIG"
  OIFS="$IFS"
  IFS=":" 
  while read line
  do
    ARR=(${line})
    # 前後の空白とタブを削除
    val=$(echo ${ARR[1]} | sed -e 's/^[ \t]*//' | sed -e 's/[ \t]*$//')
    if [[ "$line" =~ "Contamination" ]] || [[ "$line" =~ "contamination" ]]; then
      CONTAMINATION=${val} 
      echo "MAX CONTAMINATION:" $CONTAMINATION
      qc_summary_options+=("--max-contamination" "${CONTAMINATION}") 
    elif [[ "$line" =~ "Completeness" ]] || [[ "$line" =~ "completeness" ]]; then
      COMPLETENESS=${val} 
      echo "MINIMUM COMPLETENESS:" $COMPLETENESS
      qc_summary_options+=("--min-completeness" "${COMPLETENESS}")
    elif [[ "$line" =~ "Coverage" ]] || [[ "$line" =~ "coverage" ]]; then
      COVERAGE=${val} 
      echo "COVERAGE THRESHOLD:" $COVERAGE
      qc_summary_options+=("--coverage-thresh" "${COVERAGE}")
    elif [[ "$line" =~ "Contig" ]] || [[ "$line" =~ "contig" ]]; then
      CONTIGS=${val} 
      echo "MAX CONTIGS:" $CONTIGS
      qc_summary_options+=("--max-contigs" "${CONTIGS}")
    elif [[ "$line" =~ "Max" ]] || [[ "$line" =~ "max" ]]; then
      MAX_SIZE=${val} 
      echo "MAX GENOME SIZE:" $MAX_SIZE
      qc_summary_options+=("--max-gen-size" "${MAX_SIZE}")
    elif [[ "$line" =~ "Min" ]] || [[ "$line" =~ "min" ]]; then
      MIN_SIZE=${val} 
      echo "MIN GENOME SIZE:" $MIN_SIZE
      qc_summary_options+=("--min-gen-size" "${MIN_SIZE}")
    else
      SIZE=${val} 
      echo "GENOME_SIZE:" $SIZE
    fi
  done < $CONFIG
  IFS="$OIFS"
fi

# crate and move into output directory
mkdir "${OUTPUT_DIR}"
cp ${LIST} "${OUTPUT_DIR}"
cp ${DIR}/grape_version.txt "${OUTPUT_DIR}"    # Copy version information to current

cd "${OUTPUT_DIR}"

# FASTP result directory
FASTP_DIR="${OUTPUT_DIR}/fastp"
mkdir "${FASTP_DIR}"

bash ${DIR}/cov_calculator_fastp.sh "${INITIAL_DIR}" "${LIST}" "${FASTP_DIR}" "${THREAD}" "${SIZE}"


# get new list of over coverage
if [ -n "$COVERAGE" ]; then
  source activate perl-bioperl-core
  perl ${DIR}/get_cov_over40.pl "${OVER_COVERAGE_STRAIN_LIST}" "$COVERAGE"
  conda deactivate
else
  extract_strain_names "${LIST}" "${OVER_COVERAGE_STRAIN_LIST}"
fi

# Assembly and QC
num_over_coverage=$(cat "${OVER_COVERAGE_STRAIN_LIST}" | wc -l)
if [ "$num_over_coverage" -gt 0 ]; then
  assembly_and_qc
  qc_summary_options+=("--checkm-result" "${CHECKM_RESULT}")
  qc_summary_options+=("--assembly-summary" "${ASSEMBLY_SUMMARY_FILE}")
fi

source activate original
COVERAGE_FILE="${OUTPUT_DIR}/coverage.txt"
python ${DIR}/write_grape_qc.py "${COVERAGE_FILE}" ${qc_summary_options[@]}
conda deactivate

# remove fastp directory
if [ $CLEANUP = "1" ]; then
  rm -r ${FASTP_DIR}
  rm -r ${CHECKM_DIR}
fi
