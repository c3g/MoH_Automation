#!/bin/bash

set -e

# Default values
INPUT_FILE=""

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_FILE="$2"; shift ;;
        -h|--help)
            echo "Usage: $0 -i <input_json>"
            exit 0
            ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$INPUT_FILE" ]]; then
    echo "Error: Input file is required. Use -i <input_json>"
    exit 1
fi

# Source the pt-cli virtual environment from abacus
source /lb/project/mugqic/projects/MOH/project_tracking_cli/venv/bin/activate

TMP_FILE=$(mktemp)
cp "$INPUT_FILE" "$TMP_FILE"

specimen_count=$(jq '.specimen | length' "$TMP_FILE")

for ((i=0; i<specimen_count; i++)); do
    specimen_name=$(jq -r ".specimen[$i].specimen_name" "$TMP_FILE")
    specimen_id=$(pt-cli route "/project/1/specimens?json={\"specimen_name\":\"$specimen_name\"}" | jq '.DB_ACTION_OUTPUT[0].id')

    # Skip if specimen not found in DB
    if [[ -z "$specimen_id" ]]; then
        continue
    fi

    sample_ids=$(pt-cli route "/project/1/specimens/$specimen_id" | jq -r '.samples[]')

    for ((j=0; ; j++)); do
        sample_path=".specimen[$i].sample[$j]"
        sample_name=$(jq -r "$sample_path.sample_name" "$TMP_FILE" 2>/dev/null)

        [[ "$sample_name" == "null" ]] && break

        # Extract prefix and suffix type from current sample name
        prefix=$(echo "$sample_name" | sed -E 's/^(.*-)[0-9]+(DT|DN|RT)$/\1/')
        suffix=$(echo "$sample_name" | grep -oE '(DT|DN|RT)$')

        # Skip if extraction failed
        if [[ -z "$prefix" || -z "$suffix" ]]; then
            continue
        fi

        for sample_id in $sample_ids; do
            db_sample_name=$(pt-cli route "/project/1/samples/$sample_id" | jq -r '.name')

            # Skip if names are identical
            if [[ "$sample_name" == "$db_sample_name" ]]; then
                continue
            fi

            # Skip if alias already matches
            existing_alias=$(jq -r "$sample_path.sample_alias // empty" "$TMP_FILE")
            if [[ "$sample_name" == "$existing_alias" ]]; then
                continue
            fi

            # Match if db sample starts with same prefix and ends with same suffix
            if [[ "$db_sample_name" =~ ^${prefix}[0-9]+${suffix}$ ]]; then
                jq "$sample_path.sample_name = \"$db_sample_name\" | $sample_path.sample_alias = [\"$sample_name\""] "$TMP_FILE" > tmp.json && mv tmp.json "$TMP_FILE"
                echo "WARNING: Sample '$sample_name' has been renamed to match existing one '$db_sample_name'"
                break
            fi
        done
    done
done

mv "$TMP_FILE" "$INPUT_FILE"
