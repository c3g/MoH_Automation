#!/usr/bin/env bash
# fix_bam_paths.sh
# Safely rewrites BAM column (last column) without disturbing any other fields.

input="$1"

awk -F'\t' -v OFS='\t' '
NR==1 {
    # Print header unchanged
    print $0
    next
}

{
    # BAM is ALWAYS the last column
    bamCol = NF

    sample  = $1
    readset = $2

    # Readset prefix (before first dot)
    nrs = split(readset, parts, /\./)
    rs_prefix = parts[1]

    # Only change BAM when sample != prefix and Readset contains a dot
    if (nrs >= 2 && sample != rs_prefix) {

        # Extract filename from existing BAM
        nb = split($bamCol, bparts, "/")
        bamfile = bparts[nb]

        # Fallback filename if empty
        if (bamfile == "" || $bamCol ~ /\/$/)
            bamfile = rs_prefix ".sorted.bam"

        # Rewrite BAM to raw_reads/<prefix>/<filename>
        $bamCol = "/lb/project/mugqic/projects/MOH/MAIN/raw_reads/" rs_prefix "/" bamfile
    }

    print $0
}
' "$input"
