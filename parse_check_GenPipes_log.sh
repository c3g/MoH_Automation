#!/bin/bash

# Function to display usage information
usage() {
  echo "Usage: $0 -l log_file.txt -o output_path"
  exit 1
}

# Parse command-line options using getopts
while getopts ":l:o:" opt; do
  case $opt in
    l)
      LOG_FILE="$OPTARG"
      ;;
    o)
      OUTPUT_PATH="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Check if both log file and output path are provided
if [ -z "$LOG_FILE" ] || [ -z "$OUTPUT_PATH" ]; then
  usage
fi

# Read the log file line by line
while IFS= read -r line; do
  # Skip log_report.py verbosity lines
  if [[ $line == *"ERROR:__main__:"* || $line == "WARNING:" ]]; then
    continue
  # Extract the pipeline, protocol, patient, and job file using grep and awk
  elif echo "$line" | grep -q "Checking"; then
    pipeline=$(echo "$line" | awk -F'[/.]' '{print $(NF-8)}')
    if [[ $pipeline == "RnaSeqLight" ]]; then
      protocol=""
    else
      protocol=$(echo "$line" | awk -F'[/.]' '{print $(NF-7)}')
    fi
    patient=$(echo "$line" | awk -F'[/.]' '{print $(NF-6)}')

  # Handle lines with 'ERROR: GenPipes json file is empty'
  elif echo "$line" | grep -q "ERROR: GenPipes json file is empty"; then
    pipeline=$(echo "$line" | awk -F'[/.]' '{print $(NF-4)}')
    if [[ $pipeline == "RnaSeqLight" ]]; then
      protocol=""
    else
      protocol=$(echo "$line" | awk -F'[._]' '{print $(NF-4)}')
    fi
    json_path=$(echo "$line" | awk -F'Cf. ' '{print $2}')
    empty_json_file="${OUTPUT_PATH}/empty_json.${pipeline}.${protocol}.txt"

    # Append the json path to the empty_json file
    echo "${json_path}" >> "${empty_json_file}"

  # Handle lines with 'WARNING: Missing files'
  elif echo "$line" | grep -q "WARNING: Missing files"; then
    pipeline=$(echo "$line" | awk -F'[/.]' '{print $(NF-8)}')
    if [[ $pipeline == "RnaSeqLight" ]]; then
      missing_job_file="${OUTPUT_PATH}/missing_files.${pipeline}.txt"
      job_list=$(echo "$line" | awk -F'in ' '{print $2}' | awk '{print $1}' | sed 's/.$//')
      echo "${job_list}" >> "${missing_job_file}"
    else
      protocol=$(echo "$line" | awk -F'[._]' '{print $(NF-8)}')
      if [[ $pipeline == $protocol ]]; then
        pipeline=$(echo "$line" | awk -F'[/.]' '{print $(NF-9)}')
      fi
      job_list=$(echo "$line" | awk -F'in ' '{print $2}' | awk '{print $1}' | sed 's/.$//')
      missing_job_file="${OUTPUT_PATH}/missing_files.${pipeline}.${protocol}.txt"

      # Append the job list to the missing_job file
      echo "${job_list}" >> "${missing_job_file}"
    fi

  # Extract the status and job file using grep and awk
  elif echo "$line" | grep -q -E "INFO|ERROR|WARNING|SUCCESS"; then
    status=$(echo "$line" | awk '{print $1}' | tr -d ':')

    # Determine the status file name based on the status
    case $status in
      INFO)
        status_file="${OUTPUT_PATH}/running.${pipeline}.${protocol}.txt"
        ;;
      ERROR)
        status_file="${OUTPUT_PATH}/error.${pipeline}.${protocol}.txt"
        ;;
      WARNING)
        status_file="${OUTPUT_PATH}/warning.${pipeline}.${protocol}.txt"
        ;;
      SUCCESS)
        status_file="${OUTPUT_PATH}/success.${pipeline}.${protocol}.txt"
        ;;
    esac

    # Append the patient information to the appropriate file based on status
    echo "${patient}" >> "${status_file}"
  fi
done < "$LOG_FILE"
