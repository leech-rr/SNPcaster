#!/bin/bash -eu

CMDNAME=`basename $0`

function _usage() {
  cat << __EOF__

  Usage: $(basename "$0") [OPTION]...

  -i      Accession Number List (Required!)
  -t      Number of Threads     (Default: 8)
  -o      Output directory path (Default: .)

  [File Format Specifications]
  * Accession number list(-i option) is a text file containing accession numbers, one per line.
    Example:
      DRR147067
      DRR147069
      ...
__EOF__
}

ACCESSION_LIST=""
THREADS=8
OUTPUT_DIR="."
while getopts "i:t:o:" optKey; do
  case "$optKey" in
  i)
    ACCESSION_LIST=${OPTARG}
    ;;
  t)
    THREADS=${OPTARG}
    ;;
  o)
    OUTPUT_DIR=${OPTARG}
    ;;
  *)
    _usage
    exit 1
    ;;
  esac
done

# Check parameters
if [ -z "$ACCESSION_LIST" ]; then
  echo "Error: Accession number list(-i option) is required."
  _usage
  exit 1
fi
if [ ! -r "$ACCESSION_LIST" ]; then
  echo "Error: Accession number list(-i option) [$ACCESSION_LIST] does not exist or is not readable."
  _usage
  exit 1
fi

# Check if threads is a number
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
  echo "Error: Number of threads(-t option) [$THREADS] must be a positive integer."
  _usage
  exit 1
fi

# Check if output dir is a dir
if [ ! -d "$OUTPUT_DIR" ]; then
  echo "Error: Output directory(-o option) [$OUTPUT_DIR] does not exist or is not a directory."
  _usage
  exit 1
fi

# Display parameters
echo "Accession number list: $ACCESSION_LIST"
echo "Number of threads: $THREADS"
echo "Output directory: $OUTPUT_DIR"
echo

source activate sra-tools
for strain in `cat "${ACCESSION_LIST}"`; do
  echo "Fetching: $strain ..."
  fasterq-dump "${strain}" -e "${THREADS}" -3 -O "${OUTPUT_DIR}"
  echo "Done: $strain"
  echo
done
conda deactivate

# Compress the output files
cd "${OUTPUT_DIR}"
echo "Compressing output files into .gz format ..."
pigz -v *fastq


# Convert to simple format
echo "Renaming output files to simple format ..."
rename -v "s/_1.fastq.gz/_R1.fastq.gz/" *_1.fastq.gz
rename -v "s/_2.fastq.gz/_R2.fastq.gz/" *_2.fastq.gz
