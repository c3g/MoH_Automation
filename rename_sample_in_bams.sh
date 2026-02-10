#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: rename_sample_in_bams.sh -o OLD_SAMPLE -n NEW_SAMPLE -d WORKDIR [-f|--force]

Options:
  -o OLD_SAMPLE   Old sample name to replace (value of SM: in @RG lines).
  -n NEW_SAMPLE   New sample name to set (value of SM: in @RG lines).
  -d WORKDIR      Directory containing symlinked *.bam files to replace.
  -f, --force     Replace BAM even if it is not a symlink (overwrite in place).
  -h              Show this help.

Behavior:
  - For each *.bam in WORKDIR:
      * Reads from the existing file (symlink OK)
      * Creates a temporary BAM with SM:NEW in @RG lines
      * If the *.bam is a symlink, removes the symlink and installs the new real BAM at the same path
      * If the *.bam is NOT a symlink:
          - Default: left untouched (script skips deletion and refuses overwrite)
          - With --force/-f: the file is replaced in place
      * If *.bam.bai is a symlink, removes it; then creates a new real .bai via: samtools index -b
  - Replacement only touches headers (@RG SM:), alignments are not modified.

Requires:
  - samtools, awk

Example:
  rename_sample_in_bams.sh \
    -o "sample_name_2" \
    -n "sample_name_1" \
    -d /path/to/raw_reads/sample_name_2 \
    -f
USAGE
}

# --- Preprocess long options to short equivalents ---
# Supports --force as -f
args=()
for a in "$@"; do
  case "$a" in
    --force) args+=("-f") ;;
    *) args+=("$a") ;;
  esac
done
set -- "${args[@]}"

# --- Parse options ---
old=""
new=""
workdir=""
force=false
while getopts ":o:n:d:fh" opt; do
  case "$opt" in
    o) old="$OPTARG" ;;
    n) new="$OPTARG" ;;
    d) workdir="$OPTARG" ;;
    f) force=true ;;
    h) usage; exit 0 ;;
    :) echo "Error: Option -$OPTARG requires an argument." >&2; usage; exit 2 ;;
    \?) echo "Error: Invalid option -$OPTARG" >&2; usage; exit 2 ;;
  esac
done

# --- Validate inputs ---
if [[ -z "$old" || -z "$new" || -z "$workdir" ]]; then
  echo "Error: -o, -n, and -d are required." >&2
  usage
  exit 2
fi
if [[ ! -d "$workdir" ]]; then
  echo "Error: WORKDIR does not exist or is not a directory: $workdir" >&2
  exit 2
fi

# --- Check dependencies ---
for cmd in samtools awk; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Error: Required command not found in PATH: $cmd" >&2
    exit 127
  fi
done

# --- Process BAM files ---
shopt -s nullglob
bam_files=( "$workdir"/*.bam )
shopt -u nullglob

if (( ${#bam_files[@]} == 0 )); then
  echo "No BAM files found in: $workdir"
  exit 0
fi

echo "Found ${#bam_files[@]} BAM file(s) in: $workdir"
echo "Renaming SM:${old} -> SM:${new}"
$force && echo "Force mode: will overwrite non-symlink BAMs."

# Iterate and process each BAM
for bam in "${bam_files[@]}"; do
  base="$(basename "$bam")"
  bai="${bam%.bam}.bai"
  tmp="${bam%.bam}.tmp.bam"

  echo "  Processing: $base"

  # Preflight: only proceed if header actually contains @RG lines with SM:OLD
  if ! samtools view -H "$bam" | grep -qE '^@RG.*[[:space:]]SM:'"$old"'([[:space:]]|$)'; then
    echo "  No SM:${old} found in @RG header lines — skipping (no changes needed)."
    continue
  fi

  # Create a temporary BAM with header SM replaced (only within @RG lines).
  # Pattern targets '@RG' lines and the SM value: SM:OLD(<tab>|end) in SM:NEW\1
  if ! samtools view -H "$bam" \
      | awk -v old="$old" -v new="$new" 'BEGIN{FS=OFS="\t"}
        /^@RG/ {
          for (i=1; i<=NF; i++) {
            if ($i ~ /^SM:/) {
              sm = substr($i, 4)
              if (sm == old) $i = "SM:" new
            }
          }
        }
        {print}' \
      | samtools reheader - "$bam" > "$tmp"; then
    echo "Error: reheader failed for $base" >&2
    rm -f "$tmp" 2>/dev/null || true
    exit 1
  fi

  if [[ -L "$bam" ]]; then
    # Original BAM is a symlink: replace with real file
    rm -f "$bam"
    mv -f "$tmp" "$bam"
    echo "  Replaced BAM symlink with real file: $base"
  else
    if $force; then
      mv -f "$tmp" "$bam"
      echo "  Replaced non-symlink BAM in place (force): $base"
    else
      echo "  Skipped: $base is not a symlink (won't overwrite). Temporary output at: $tmp"
      continue
    fi
  fi

  # If BAI symlink exists, remove it before indexing to avoid writing through the symlink.
  if [[ -L "$bai" ]]; then
    rm -f "$bai"
    echo "  Removed BAI symlink: $(basename "$bai")"
  fi

  # Create new BAI index (real file) for the new BAM.
  if ! samtools index -b "$bam" "$bai"; then
    echo "Error: indexing failed for $base" >&2
    exit 1
  fi
  echo "  Wrote new index: $(basename "$bai")"
done

echo "All done."
