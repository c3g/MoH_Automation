#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: rename_sample_in_bams.sh -o OLD_SAMPLE -n NEW_SAMPLE -d WORKDIR

Options:
  -o OLD_SAMPLE   Old sample name to replace (value of SM: in @RG lines).
  -n NEW_SAMPLE   New sample name to set (value of SM: in @RG lines).
  -d WORKDIR      Directory containing symlinked *.bam files to replace.
  -h              Show this help.

Behavior:
  - For each *.bam in WORKDIR:
      * Reads from the existing file (symlink OK)
      * Creates a temporary BAM with SM:NEW in @RG lines
      * If the *.bam is a symlink, removes the symlink and installs the new real BAM at the same path
      * If *.bam.bai is a symlink, removes it and creates a new real .bai via samtools index -b
  - If a *.bam is NOT a symlink, it is left untouched (script skips deletion and refuse overwrite).

Requires:
  - samtools, sed

Example:
  rename_sample_in_bams.sh \
    -o "sample_name_2" \
    -n "sample_name_1" \
    -d /path/to/raw_reads/sample_name_2
USAGE
}

# --- Parse options ---
old=""
new=""
workdir=""
while getopts ":o:n:d:h" opt; do
  case "$opt" in
    o) old="$OPTARG" ;;
    n) new="$OPTARG" ;;
    d) workdir="$OPTARG" ;;
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
for cmd in samtools sed; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Error: Required command not found in PATH: $cmd" >&2
    exit 127
  fi
done

# --- Process BAM files safely (handles spaces with NUL separation) ---
shopt -s nullglob
bam_files=( "$workdir"/*.bam )
shopt -u nullglob

if (( ${#bam_files[@]} == 0 )); then
  echo "No BAM files found in: $workdir"
  exit 0
fi

echo "Found ${#bam_files[@]} BAM file(s) in: $workdir"
echo "Renaming SM:${old} -> SM:${new}"

# Iterate and process each BAM
for bam in "${bam_files[@]}"; do
  base="$(basename "$bam")"
  bai="${bam%.bam}.bai"
  tmp="${bam%.bam}.tmp.bam"

  echo "  Processing: $base"

  # Create a temporary BAM with header SM replaced (only within @RG lines).
  # Pattern targets '@RG' lines and the SM value: SM:OLD(<tab>|end) in SM:NEW\1
  if ! samtools view -H "$bam" \
      | sed -E "s/^(@RG[^\t]*\tSM:)${old}(\t|$)/\1${new}\2/g" \
      | samtools reheader - "$bam" > "$tmp"; then
    echo "Error: reheader failed for $base" >&2
    rm -f "$tmp" 2>/dev/null || true
    exit 1
  fi

  # Ensure the original BAM is a symlink before replacing.
  if [[ -L "$bam" ]]; then
    # Remove BAM symlink and move temp file.
    rm -f "$bam"
    mv -f "$tmp" "$bam"
    echo "  Replaced BAM symlink with real file: $base"
  else
    echo "  Skipped: $base is not a symlink (won't overwrite). Temporary output at: $tmp"
    # Optionally continue; or you could bail out. Here we continue and leave tmp for manual check.
    continue
  fi

  # If BAI symlink exists, remove it before indexing to avoid writing through the symlink.
  if [[ -L "$bai" ]]; then
    rm -f "$bai"
    echo "  Removed BAI symlink: $(basename "$bai")"
  fi

  # Create new BAI index (real file) for the new BAM.
  if ! samtools index -b "$bam"; then
    echo "Error: indexing failed for $base" >&2
    exit 1
  fi
  echo "  Wrote new index: $(basename "$bai")"

done

echo "All done."
