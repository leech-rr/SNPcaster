#!/usr/bin/env bash
#
# objective
#   Batch-run amrfinder for a list of strains with flexible options.
# usage
#   bash amrfinder_batch.sh -i strains.txt -e .fna -t 8 --plus -s "Escherichia coli"
# dependency
#   amrfinder (NCBI AMRFinderPlus), GNU coreutils, bash >= 4
#
# Notes
#   * Output folder: amrfinder_YYYYMMDD_HHMMSS/
#   * Logs are collected under amrfinder_YYYYMMDD_HHMMSS/log/
#   * Output filenames are: S_plus_out (if --plus) or S_out (if not).
#
set -euo pipefail

show_help() {
  cat <<'EOF'
Usage: bash amrfinder_batch.sh -i LIST -e EXT [--plus] [-t THREADS] [-s SPECIES]

  -i, --input         Text file of strain names (one per line) [required]
  -e, --ext           Genome file extension (e.g., .fna or fna) [required]
  -t, --threads       Number of threads for amrfinder (optional, default: 8)
      --plus          Add --plus option to amrfinder (optional)
  -s, --species       Species string for amrfinder -O (e.g., "Escherichia") (optional)
                      Use --list-species to see available organisms
  -h, --help          Show this help
      --list-species  Show list of available species/organisms for -s option
EOF
}

list_species() {
  source activate amrfinder
  if command -v amrfinder &>/dev/null; then
    amrfinder -l 2>/dev/null | sed '/^$/d; s/--organism/-s\/--species/g; s/^/[INFO] /' || echo "[WARN] Could not retrieve species list from amrfinder"
  else
    echo "[ERROR] amrfinder not found in PATH. Please activate amrfinder environment first."
    conda deactivate
    exit 1
  fi
  conda deactivate
}

require_argument() {
  local option="$1"
  local value="$2"
  if [[ -z "$value" || "$value" =~ ^- ]]; then
    echo "[ERROR] $option requires an argument" >&2
    show_help; exit 1
  fi
}

# Defaults
INPUT_LIST=""
EXT=""
THREADS="8"
PLUS="false"
SPECIES=""

if [[ $# -eq 0 ]]; then
  show_help
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      require_argument "-i/--input" "${2:-}"
      INPUT_LIST="$2"; shift 2;;
    -e|--ext)
      require_argument "-e/--ext" "${2:-}"
      EXT="$2"; shift 2;;
    -t|--threads)
      require_argument "-t/--threads" "${2:-}"
      THREADS="$2"; shift 2;;
    --plus)
      PLUS="true"; shift;;
    -s|--species)
      require_argument "-s/--species" "${2:-}"
      SPECIES="$2"; shift 2;;
    -h|--help)
      show_help; exit 0;;
    --list-species)
      list_species
      exit 0;;
    --) shift; break;;
    -*)
      echo "[ERROR] Unknown option: $1" >&2
      show_help; exit 1;;
    *)
      echo "[ERROR] Unexpected positional argument: $1" >&2
      show_help; exit 1;;
  esac
done

if [[ -z "$INPUT_LIST" || -z "$EXT" ]]; then
  echo "[ERROR] -i and -e are required." >&2
  show_help
  exit 1
fi

if [[ ! -f "$INPUT_LIST" ]]; then
  echo "[ERROR] Input list not found: $INPUT_LIST" >&2
  exit 1
fi

# Normalize extension
if [[ "$EXT" != .* ]]; then
  EXT=".$EXT"
fi

# Create timestamped folders
TS_DIR="amrfinder_$(date +%Y%m%d_%H%M%S)"
LOG_DIR="$TS_DIR/log"
mkdir -p "$LOG_DIR"
echo "[INFO] Created output folder: $TS_DIR"
echo "[INFO] Created log folder: $LOG_DIR"

cleanup() {
  shopt -s nullglob
  for f in *_out *_plus_out; do
    [[ -e "$f" ]] || continue
    mv -f "$f" "$TS_DIR/" 2>/dev/null || true
  done
  for f in *.log; do
    [[ -e "$f" ]] || continue
    mv -f "$f" "$LOG_DIR/" 2>/dev/null || true
  done
}
trap cleanup EXIT

# Validation
source activate amrfinder
if ! command -v amrfinder &>/dev/null; then
  echo "[ERROR] amrfinder not found in PATH." >&2
  exit 1
fi

# Process each strain
while IFS=$'\n' read -r strain || [[ -n "$strain" ]]; do
  [[ -z "$strain" || "$strain" =~ ^# ]] && continue

  genome="${strain}${EXT}"

  if [[ ! -f "$genome" ]]; then
    echo "[WARN] Genome file not found for strain '${strain}': $genome" >&2
    continue
  fi

  log="$LOG_DIR/${strain}.log"

  # Build amrfinder command arguments
  amrfinder_args=(
    "-n" "$genome"
    "--threads" "$THREADS"
  )
  if [[ -n "$SPECIES" ]]; then
    amrfinder_args+=("--organism" "$SPECIES")
  fi
  out="${strain}_out.tsv"
  if [[ "$PLUS" == "true" ]]; then
    amrfinder_args+=("--plus")
    out="${strain}_plus_out.tsv"
  fi

  echo "[INFO] Running amrfinder on $genome -> $out"
  amrfinder "${amrfinder_args[@]}" > "$out" 2> "$log"

  mv -f "$out" "$TS_DIR/"
done < "$INPUT_LIST"

echo "[INFO] Done. Results in: $TS_DIR ; logs in: $LOG_DIR"

conda deactivate
