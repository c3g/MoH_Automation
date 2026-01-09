#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: rename_sample_in_bams.sh -o OLD_SAMPLE -n NEW_SAMPLE -d WORKDIR

Options:
  -o OLD_SAMPLE   Old sample name to replace (value of SM: in @RG lines).
  -n NEW_SAMPLE   New sample name to set (value of SM: in @RG lines).
  -d WORKDIR      Directory containing *.bam files to process.
  -h              Show this help.

Example:
  rename_sample_in_bams.sh \
    -o "MoHQ-JG-26-1-MD180207-2DN" \
    -n "MoHQ-JG-26-1-MD180207-1DN" \
    -d /data/project/run42

Notes:
  - Output files are named: <basename>.renamed.bam and <basename>.renamed.bam.bai
  - Exact match replacement is enforced for SM:OLD_SAMPLE only.
  - Requires: samtools, sed
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
mapfile -d '' bam_files < <(find "$workdir" -maxdepth 1 -type f -name '*.bam' -print0)
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
  out="${bam%.bam}.renamed.bam"

  echo "Processing: $base"
  # Reheader with exact SM match; global replacement across header
  # \b ensures we hit the exact token right after 'SM:' and not substrings.
  # Using sed -E for extended regex and escaping backslashes properly.
  if ! samtools view -H "$bam" \
      | sed -E "s/\bSM:${old}\b/SM:${new}/g" \
      | samtools reheader - "$bam" > "$out"; then
    echo "Error: reheader failed for $base" >&2
    exit 1
  fi

  # Build BAI index for the renamed BAM
  if ! samtools index -b "$out"; then
    echo "Error: indexing failed for $(basename "$out")" >&2
    exit 1
  fi

  # Quick verification (optional): show updated RG SM values
  # samtools view -H "$out" | grep '^@RG' || true

  echo "Wrote: $(basename "$out") and $(basename "$out").bai"
done

echo "All done."
