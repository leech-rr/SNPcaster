#!/bin/bash

# This script only has functions for SNPcaster.
# The main script is "snpcaster.sh".

SNPCASTER_SRC_DIR="${BAC1}"/snpcaster
CORE_RCM_SRC_DIR="${SNPCASTER_SRC_DIR}"/core_rcm
source "${BAC1}"/grape_qc_assembly/strain_list_functions.sh

# Change the following variables if you want to change the input/output directory names.
BACTSNP_OUTDIR=1_results_bactsnp
SNIPPY_OUTDIR=2_snippy_results
RESULTS_WITHOUT_GUBBINS_OUTDIR=3_results_without_gubbins
GUBBINS_OUTDIR=4_results_gubbins
RESULTS_WITH_GUBBINS_OUTDIR=5_results_with_gubbins
EXT=.fasta
TITLE4DIST="strain1\tstrain2\tpairwise_snp"
FINAL_CORE_REGION_FILE_NAME=core_region_final.tsv

function create_pairsnp_matrix() {
  local final_snp=$1
  local output_path=$2
  local THREAD=$3
  local temp_path="${output_path}_temp"

  source activate pairsnp
  pairsnp -t "${THREAD}" "${final_snp}" > "${temp_path}"
  conda deactivate
  # create a header from column 1
  awk 'BEGIN {ORS=""} {print "\t" $1 } END {print "\n"}' "${temp_path}" > "${output_path}"
  # add a header into output
  cat "${temp_path}" >> "${output_path}"
  rm "${temp_path}"
}

function rename_fasta_extension() {
  local TARGET_DIR="${1%/}/"
  find "${TARGET_DIR}" -type f -name "*.fa" -exec rename -v 's/\.fa$/'"${EXT}"'/' {} +
  find "${TARGET_DIR}" -type f -name "*.fas" -exec rename -v 's/\.fas$/'"${EXT}"'/' {} +
}

function get_version() {
  cat "${SNPCASTER_SRC_DIR}"/snpcaster_version.txt
}

function count_lines_without_header() {
  local INPUT="$1"
  local header_lines=1
  if [ -n "$2" ]; then
    header_lines="$2"
  fi
  tail -n +$((header_lines + 1)) "${INPUT}" | wc -l
}

function count_core_genome_length() {
  local CORE_GENOME_TSV="$1"
  if [ ! -r "${CORE_GENOME_TSV}" ]; then
    echo "0"
    return
  fi
  source activate original
  python "${CORE_RCM_SRC_DIR}"/core_genome_calc.py "${CORE_GENOME_TSV}" # output: core_genome.txt
  conda deactivate
  local core_genome_count=0
  if [ -r "core_genome.txt" ]; then
    core_genome_count=$(cat core_genome.txt)
    rm core_genome.txt
  fi
  echo "${core_genome_count}"
}

function prepare4snpcaster() {
  local NOT_EXIST_LIST=missing_list
  local BACTSNP_LIST=$1
  local FASTQ_LIST=$2
  local LIST_FULL=$3
  local OUTPUT_DIR=$4
  local BACTSNP_SKIP_LIST=$5
  local REFERENCE_FILE=$6
  local MASK_FILE=$7

  [ -f "${BACTSNP_LIST}" ] && rm -f "${BACTSNP_LIST}"
  [ -f "${BACTSNP_SKIP_LIST}" ] && rm -f "${BACTSNP_SKIP_LIST}"
  [ -f "${NOT_EXIST_LIST}" ] && rm -f "${NOT_EXIST_LIST}"

  touch $BACTSNP_SKIP_LIST
  local strain
  for strain in $(cat "${LIST_FULL}"); do
    # Skip if pre-executed SNP folder exists
    if [ -d "${strain}" ]; then
      cp -pvr "${strain}" "${OUTPUT_DIR}"/
      rename 's/\.fasta$/.fa/' "${OUTPUT_DIR}/${strain}"/*.fasta
      echo "${strain}" >> "${BACTSNP_SKIP_LIST}"
      continue
    fi

    # User must specify fastq list to execute BactSNP
    if [ -z "${FASTQ_LIST}" ]; then
      echo "Error: Please specify fastq list(-f option) for ${strain} to execute BactSNP." >&2
      _usage
      exit 1
    fi

    # Find fastq file pair for the strain
    local fastq_line=$(grep -P "^${strain}\t" "${FASTQ_LIST}" | tr -d '\n')

    if [ -z "${fastq_line}" ]; then
      echo "${strain}" >> "${NOT_EXIST_LIST}"
      echo "Error: Could not find ${strain} in [${FASTQ_LIST}]." >&2
      continue
    fi

    IFS=$'\t' read -r -a columns <<< "${fastq_line}"
    local read1_file_name=${columns[1]}
    local read2_file_name=${columns[2]}

    if [ ! -r "${read1_file_name}" ]; then
      echo "${strain}" >> "${NOT_EXIST_LIST}"
      echo "Error: Could not find or read ${read1_file_name} file (strain=${strain}) written in $FASTQ_LIST." >&2
      continue
    fi
    if [ ! -r "${read2_file_name}" ]; then
      echo "${strain}" >> "${NOT_EXIST_LIST}"
      echo "Error: Could not find or read $read2_file_name file (strain=${strain}) written in $FASTQ_LIST." >&2
      continue
    fi

    # BactSNP will be executed in OUTPUT_DIR directory,
    # so ../ is needed for file names
    echo -e "${strain}""\t../""${read1_file_name}""\t../""${read2_file_name}" >> "${BACTSNP_LIST}"
  done

  if [ -f "${NOT_EXIST_LIST}" ]; then
    echo ""
    echo "*********** The fastq file(s) does not exist. See missing_list. **************"
    cat "${NOT_EXIST_LIST}"
    mv "${NOT_EXIST_LIST}" "${OUTPUT_DIR}"/
    #すでにコピー済みの株名フォルダを削除
    for strain in $(
      cat "${BACTSNP_SKIP_LIST}"
    ); do
      rm -r "${OUTPUT_DIR}"/"${strain}"/
    done
    exit 1
  fi

  cp -v "${LIST_FULL}" "${OUTPUT_DIR}"/ || exit 1
  cp -v "${REFERENCE_FILE}" "${OUTPUT_DIR}"/ || exit 1
  if [ -n "${MASK_FILE}" ]; then
    cp -v "${MASK_FILE}" "${OUTPUT_DIR}"/ || exit 1
  fi
  [ -f "${BACTSNP_LIST}" ] && mv -v "${BACTSNP_LIST}" "${OUTPUT_DIR}"/
  [ -f "${BACTSNP_SKIP_LIST}" ] && mv -v "${BACTSNP_SKIP_LIST}" "${OUTPUT_DIR}"/
}

function process_bactsnp() {
  local BACTSNP_LIST=$1
  local REFERENCE_FILE=$2
  local ALLELE_FREQ=$3
  local THREAD=$4
  local JOBS=$5

  echo '===================BactSNP==================='
  if [ -f "${BACTSNP_LIST}" ]; then
    [ -d "${BACTSNP_OUTDIR}" ] && rm -r "${BACTSNP_OUTDIR}"
    which bactsnp
    echo "bactsnp -r ${REFERENCE_FILE} --fastq_list ${BACTSNP_LIST} --allele_freq ${ALLELE_FREQ} -t ${THREAD} -j ${JOBS} -o ${BACTSNP_OUTDIR}"
    source activate bactsnp-dependencies
    bactsnp -r "${REFERENCE_FILE}" --fastq_list "${BACTSNP_LIST}" --allele_freq "${ALLELE_FREQ}" -t "${THREAD}" -j "${JOBS}" -o "${BACTSNP_OUTDIR}"
    conda deactivate
  fi

  local PSEUDO_GENOME_DIR="${BACTSNP_OUTDIR}"/pseudo_genome
  local TITLE="#CHROM\tPOS\tID\tREF\tALT"
  source activate perl-bioperl-core
  for strain in $(extract_strain_names "${BACTSNP_LIST}"); do
    file=${PSEUDO_GENOME_DIR}/${strain}.fa
    perl "${SNPCASTER_SRC_DIR}"/bactsnp_vcf_script.pl "${REFERENCE_FILE}" "${file}"
    # Add title to the result file
    sed -i "1i ${TITLE}" result_vcf.txt
    # Create SNP folder for each strain
    mkdir -v "${strain}"
    mv -v result_vcf.txt "${strain}"/snps.vcf
    cp -v "${file}" "${strain}"/snps.aligned.fa
  done
  conda deactivate
}

function process_snippy() {
  local STRAIN_LIST=$1
  local REFERENCE_FILE=$2

  echo '===================Snippy==================='
  local my_snps=""
  for strain in $(cat "${STRAIN_LIST}"); do
    my_snps="${my_snps}"" ${strain}"
  done
  source activate snippy
  # shellcheck disable=SC2086
  snippy-core --prefix core ${my_snps} --ref "${REFERENCE_FILE}"
  conda deactivate
  mkdir "${SNIPPY_OUTDIR}"
  mv core.aln core.full.aln core.tab core.txt core.vcf "${SNIPPY_OUTDIR}"/
}

function postprocess_snippy() {
  local SNP_POSITION_FILE=snp_position.csv
  local GAP=$1
  local THREAD=$2
  local MASK_FILE=$3

  local using_mask_file=0
  if [ -n "${MASK_FILE}" ]; then
    using_mask_file=1
  fi

  echo '===================Snippy results conversion ==================='
  # substitution (core.tab --> snp_position.csv)
  sed -e "s/\\t/,/g" "${SNIPPY_OUTDIR}"/core.tab > "${SNP_POSITION_FILE}"
  sed -i -e "s/,LOCUS_TAG,GENE,PRODUCT,EFFECT//" -e "s/,,,,\$//g" "${SNP_POSITION_FILE}"

  source activate perl-bioperl-core
  local snp_conversion_options=("${SNP_POSITION_FILE}" "-c" "${GAP}")
  if [ "${using_mask_file}" -eq 1 ]; then
    snp_conversion_options+=("-d" "${MASK_FILE}")
  fi
  perl "${SNPCASTER_SRC_DIR}"/snp_conversion.pl ${snp_conversion_options[@]}
  conda deactivate

  source activate pairsnp
  pairsnp -t "${THREAD}" -s final_snp.fas > dist_final_snp.tsv
  conda deactivate
  create_pairsnp_matrix final_snp.fas dist_final_snp_matrix.tsv "${THREAD}"

  # 1行目にタイトルを追記
  sed -i "1i ${TITLE4DIST}" dist_final_snp.tsv

  # final_snpのRef除外版とnexus変換版を作成**********************************
  bash "${SNPCASTER_SRC_DIR}"/seqmagick_conversion.sh

  # Calculate core genome size
  source activate perl-bioperl-core
  perl "${CORE_RCM_SRC_DIR}"/core_extract.pl "${SNIPPY_OUTDIR}"/core.full.aln
  local final_core_region_src_file=core_region.tsv
  if [ "${using_mask_file}" -eq 1 ]; then
    perl "${CORE_RCM_SRC_DIR}"/core_remove_rep_rcm.pl core_region.tsv "${MASK_FILE}"
    final_core_region_src_file="core_region_after_masking.tsv"
    mv -v core_region_without_recombination.tsv "${final_core_region_src_file}"
  fi
  cp -pv "${final_core_region_src_file}" "${FINAL_CORE_REGION_FILE_NAME}"
  conda deactivate

  # Move core files to 3.folder.
  mkdir "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"
  mv -v "${SNP_POSITION_FILE}" snp_position_sample_only.csv snp_position_final.csv core.full.fas final_snp.fas final_snp.nex final_snp_woRef.fas final_snp_woRef.nex dist_final_snp.tsv dist_final_snp_matrix.tsv core_region.tsv "${FINAL_CORE_REGION_FILE_NAME}" "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/
  if [ "${using_mask_file}" -eq 1 ]; then
    mv -v snp_position_after_masking_sample_only.csv masked_region.csv snp_position_after_masking.csv core_region_after_masking.tsv "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/
  fi
  if [ "${GAP}" -gt 0 ]; then
    mv -v snp_position_without_clusterSNP.csv snp_position_without_clusterSNP_sample_only.csv removed_clusterSNP.csv "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/
  fi
}

function process_gubbins() {
  local THREAD=$1
  echo '========================= gubbins ==========================='
  source activate gubbins
  run_gubbins.py --threads "${THREAD}" "${SNIPPY_OUTDIR}"/core.full.aln
  conda deactivate

  mkdir "${GUBBINS_OUTDIR}"
  mv -v core.full.branch_base_reconstruction.embl \
    core.full.filtered_polymorphic_sites.fasta core.full.filtered_polymorphic_sites.phylip \
    core.full.final_tree.tre core.full.node_labelled.final_tree.tre \
    core.full.recombination_predictions.embl core.full.recombination_predictions.gff \
    core.full.summary_of_snp_distribution.vcf core.full.per_branch_statistics.csv "${GUBBINS_OUTDIR}"/
}

function postprocess_gubbins() {
  local THREAD=$1
  echo "========================= gubbins results conversion =========================="
  # create recombination.tsv
  local TMP_FILE1=extracted_position_fr_gff.txt
  local TMP_FILE2=extracted_position_fr_gff_without_overlap.txt
  local RECOMBINATION_TSV_FILE=recombination.tsv
  local SNP_POSITION_FINAL_FILE="${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/snp_position_final.csv
  local REF_ID=$(awk -v FPAT='[^,]*|"[^"]*"' 'NR==2 {print $1; exit}' "${SNP_POSITION_FINAL_FILE}")

  sed -e '1,2d' "${GUBBINS_OUTDIR}"/core.full.recombination_predictions.gff > "${TMP_FILE1}"
  awk -i inplace -F'\t' '{print $4 "\t" $5}' "${TMP_FILE1}"
  source activate perl-bioperl-core
  perl "${SNPCASTER_SRC_DIR}"/delete_overlapped_data.pl "${TMP_FILE1}"
  # remove a header
  sed -i -e '1d' "${TMP_FILE2}"
  # add REF_ID into the first column
  sed "s/^/${REF_ID}\t/" "${TMP_FILE2}" > "${RECOMBINATION_TSV_FILE}"
  # add a header
  sed -i "1s/^/Tag\tStart\tEnd\n/" "${RECOMBINATION_TSV_FILE}"
  # remove temp files
  rm "${TMP_FILE1}" "${TMP_FILE2}"
  perl "${SNPCASTER_SRC_DIR}"/snp_conversion.pl "${SNP_POSITION_FINAL_FILE}" -d "${RECOMBINATION_TSV_FILE}"
  rm -v 'snp_position_sample_only.csv' 'snp_position_final.csv'
  mv -v 'snp_position_after_masking.csv' 'snp_position_after_gubbins.csv'
  mv -v 'snp_position_after_masking_sample_only.csv' 'snp_position_after_gubbins_sample_only.csv'
  mv -v 'masked_region.csv' 'recombination_region.csv'
  mv -v 'final_snp.fas' 'final_snp_after_gubbins.fas'
  conda deactivate
  # create gubbins final snp without reference
  cp -v final_snp_after_gubbins.fas final_snp_after_gubbins_woRef.fas
  sed -i -e '1,2d' final_snp_after_gubbins_woRef.fas

  source activate pairsnp
  pairsnp -t "${THREAD}" -s final_snp_after_gubbins.fas > dist_final_snp_without_recombination.tsv
  conda deactivate
  sed -i "1i ${TITLE4DIST}" dist_final_snp_without_recombination.tsv
  create_pairsnp_matrix final_snp_after_gubbins.fas dist_final_snp_matrix_without_recombination.tsv

  # final_snp_after_gubbins.fas のRef除外版とnexus変換版を作成**********************************
  bash "${SNPCASTER_SRC_DIR}"/seqmagick_conversion.sh final_snp_after_gubbins.fas
  # calculate core genome size
  perl "${CORE_RCM_SRC_DIR}"/core_remove_rep_rcm.pl "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/"${FINAL_CORE_REGION_FILE_NAME}" "${RECOMBINATION_TSV_FILE}"
  mv -v "${FINAL_CORE_REGION_FILE_NAME%.*}"_without_recombination.tsv core_summary_after_gubbins.tsv

  mkdir -v "${RESULTS_WITH_GUBBINS_OUTDIR}"
  mv -v final_snp_after_gubbins.fas final_snp_after_gubbins.nex final_snp_after_gubbins_woRef.fas final_snp_after_gubbins_woRef.nex recombination_region.csv snp_position_after_gubbins.csv dist_final_snp_without_recombination.tsv dist_final_snp_matrix_without_recombination.tsv core_summary_after_gubbins.tsv snp_position_after_gubbins_sample_only.csv "${RESULTS_WITH_GUBBINS_OUTDIR}"/
}

function rename_fasta_extension_in_directories() {
  local target_dir=""
  for target_dir in "$@"; do
    rename_fasta_extension "${target_dir}"
  done
}

function append_report_information() {
  local REPORT=$1
  local WRITE_NO_MASKING=$2
  local WRITE_CLUSTER_SNP_REMOVAL=$3
  local WRITE_AFTER_MASKING=$4
  local WRITE_AFTER_GUBBINS=$5

  echo "" >> "${REPORT}"
  echo "Core genome size (bp) " >> "${REPORT}"
  echo "Before Gubbins : $(count_core_genome_length "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/"${FINAL_CORE_REGION_FILE_NAME}")" >> "${REPORT}"
  if [ "${WRITE_AFTER_GUBBINS}" -eq 1 ]; then
    echo "After Gubbins : $(count_core_genome_length "${RESULTS_WITH_GUBBINS_OUTDIR}"/core_summary_after_gubbins.tsv)" >> "${REPORT}"
  fi

  echo "" >> "${REPORT}"
  echo "Number of informative SNP sites (with reference) " >> "${REPORT}"
  if [ "$WRITE_NO_MASKING" -gt 0 ]; then
    echo "No masking: $(count_lines_without_header "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/snp_position.csv)" >> "${REPORT}"
  fi
  if [ "$WRITE_CLUSTER_SNP_REMOVAL" -gt 0 ]; then
    echo "After cluster SNP removal: $(count_lines_without_header "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/snp_position_without_clusterSNP.csv)" >> "${REPORT}"
  fi
  if [ "$WRITE_AFTER_MASKING" -gt 0 ]; then
    echo "After masking: $(count_lines_without_header "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/snp_position_after_masking.csv)" >> "${REPORT}"
  fi
  if [ "$WRITE_AFTER_GUBBINS" -gt 0 ]; then
    echo "After Gubbins : $(count_lines_without_header "${RESULTS_WITH_GUBBINS_OUTDIR}"/snp_position_after_gubbins.csv)" >> "${REPORT}"
  fi

  echo "" >> "${REPORT}"
  echo "Number of informative SNP sites (without reference) " >> "${REPORT}"
  if [ "$WRITE_NO_MASKING" -gt 0 ]; then
    echo "No masking: $(count_lines_without_header "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/snp_position_sample_only.csv)" >> "${REPORT}"
  fi
  if [ "$WRITE_CLUSTER_SNP_REMOVAL" -gt 0 ]; then
    echo "After cluster SNP removal: $(count_lines_without_header "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/snp_position_without_clusterSNP_sample_only.csv)" >> "${REPORT}"
  fi
  if [ "$WRITE_AFTER_MASKING" -gt 0 ]; then
    echo "After masking: $(count_lines_without_header "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/snp_position_after_masking_sample_only.csv)" >> "${REPORT}"
  fi
  if [ "$WRITE_AFTER_GUBBINS" -gt 0 ]; then
    echo "After Gubbins : $(count_lines_without_header "${RESULTS_WITH_GUBBINS_OUTDIR}"/snp_position_after_gubbins_sample_only.csv)" >> "${REPORT}"
  fi

  if [ "$WRITE_CLUSTER_SNP_REMOVAL" -gt 0 ] || [ "$WRITE_AFTER_MASKING" -gt 0 ] || [ "$WRITE_AFTER_GUBBINS" -gt 0 ]; then
    echo "" >> "${REPORT}"
    echo "Removed SNPs" >> "${REPORT}"
    if [ "$WRITE_CLUSTER_SNP_REMOVAL" -gt 0 ]; then
      echo "Cluster SNP: $(count_lines_without_header "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/removed_clusterSNP.csv)" >> "${REPORT}"
    fi
    if [ "$WRITE_AFTER_MASKING" -gt 0 ]; then
      echo "Masked region: $(count_lines_without_header "${RESULTS_WITHOUT_GUBBINS_OUTDIR}"/masked_region.csv)" >> "${REPORT}"
    fi
    if [ "$WRITE_AFTER_GUBBINS" -gt 0 ]; then
      echo "Gubbins: $(count_lines_without_header "${RESULTS_WITH_GUBBINS_OUTDIR}"/recombination_region.csv)" >> "${REPORT}"
    fi
  fi
}
