#!/usr/bin/env bash
# fix_bam_paths.sh
# Rewrites the BAM column (assumed LAST column) when Sample != prefix(Readset).
# New BAM: /lb/project/mugqic/projects/MOH/MAIN/raw_reads/<prefix>/<original_bam_filename>
#
# Usage:
#   chmod +x fix_bam_paths.sh
#   ./fix_bam_paths.sh input.tsv > output.tsv

input=$1

awk -v OFS='\t' '
BEGIN {
  RAW_BASE = "/lb/project/mugqic/projects/MOH/MAIN/raw_reads";
}
NR==1 {
  # Print header exactly as-is (TSV preserved)
  print $0;
  next;
}
{
  # BAM is the LAST column
  bamCol = NF;

  # Extract key fields by name positions based on the fixed order in your TSV:
  # Sample(1) Readset(2) LibraryType(3) RunType(4) Run(5) Lane(6) Adapter1(7) Adapter2(8) QualityOffset(9) BED(10) FASTQ1(11) FASTQ2(12) BAM(13)
  # But we only need Sample and Readset and BAM
  sample  = $1;
  readset = $2;

  # Get <prefix> from Readset (substring before the first dot)
  # If there is no dot, parts[1] will be the whole Readset
  nrs = split(readset, parts, /\./);
  rs_prefix = parts[1];

  # Only modify when Sample != rs_prefix AND we actually had a dot (nrs >= 2)
  if (rs_prefix != "" && sample != rs_prefix && nrs >= 2) {

    # Preserve original BAM filename (basename of the path)
    nb = split($bamCol, bparts, "/");
    bamfile = bparts[nb];

    # Fallback in case BAM column is empty or ends with a slash
    if (bamfile == "" || $bamCol ~ /\/$/) {
      bamfile = rs_prefix ".sorted.bam";
    }

    # Build new BAM path under raw_reads/<prefix>
    $bamCol = RAW_BASE "/" rs_prefix "/" bamfile;
  }

  print $0;
}
' "$input"