#!/usr/bin/env bash

THIS_SCRIPT=$(basename "$0")

usage() {
  echo "script usage: $THIS_SCRIPT -h [-r run_processing_json] [-d destination]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -r <run_processing_json>         Run Processing json."
  echo " -d <destination>                 Destination for the transfer (either Rorqual or Narval or Cardinal or Abacus)."
  exit 1
  }

while getopts 'hr:d:' OPTION; do
  case "$OPTION" in
    r)
      run_processing_json="$OPTARG"
      ;;
    d)
      destination="$OPTARG"
      ;;
    h)
      usage
      ;;
    ?)
      usage
      ;;
  esac
done

# mandatory arguments
if [ ! "$run_processing_json" ] || [ ! "$destination" ]; then
  echo -e "ERROR: Missing mandatory argument -r and/or -d.\n"
  usage
fi

if ! [[ $destination =~ Cardinal|Rorqual|Abacus|Narval ]]; then
    echo -e "ERROR: Invalid destination: '$destination'. It has to be either Rorqual, Narval, Cardinal or Abacus.\n"
    usage
fi

ROBOT_EP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/globus_collections.json"

destination_lowercase="${destination,,}"
DEST_EP=$(jq -r --arg dest "$destination_lowercase" '.robot_endpoints[$dest].uuid' "$ROBOT_EP")
DEST_BASE_PATH=$(jq -r --arg dest "$destination_lowercase" '.robot_endpoints[$dest].base_path[0]' "$ROBOT_EP")

if [[ $destination = Rorqual ]]; then
  # Rorqual Endpoint
  # DEST_EP='f19f13f5-5553-40e3-ba30-6c151b9d35d4'
  # Rorqual_MoH_Robot Endpoint
  # DEST_EP='b8ae5665-590d-45b9-ae83-d1dfde08e7d0'
  # Rorqual base path
  # DEST_BASE_PATH="/project/6007512/C3G/projects/MOH_PROCESSING/MAIN"
  # Rorqual main folder location
  DEST_LOC="${DEST_BASE_PATH}/MAIN/raw_reads"
  # Rorqual log file location
  DEST_LOG_LOC="${DEST_BASE_PATH}/DATABASE/log_files/transfer"
elif [[ $destination = Narval ]]; then
  # Narval main folder location
  DEST_LOC="${DEST_BASE_PATH}/MAIN/raw_reads"
  # Narval log file location
  DEST_LOG_LOC="${DEST_BASE_PATH}/DATABASE/log_files/transfer"
elif [[ $destination = Cardinal ]]; then
  # Cardinal Endpoint
  # DEST_EP='a6df16bc-4e7f-4784-9afa-8ceb7b20b7c0'
  # Cardinal_MoH_Robot Endpoint
  # DEST_EP='995863d1-97fe-4af2-bc48-c9bf91c33f08'
  # Cardinal base path
  # DEST_BASE_PATH="/project/60007/MOH"
  # Cardinal main folder location
  DEST_LOC="${DEST_BASE_PATH}/MAIN/raw_reads"
  # Cardinal log file location
  DEST_LOG_LOC="${DEST_BASE_PATH}/log_files/transfer"
elif [[ $destination = Abacus ]]; then
  # Abacus Endpoint
  DEST_EP=''
  # Abacus main folder location
  DEST_LOC="/lb/project/mugqic/projects/MOH/MAIN/raw_reads"
  # Abacus log file location
  DEST_LOG_LOC="/lb/project/mugqic/projects/MOH/log_files/transfer"
fi

# Abacus Endpoint
# ABA_EP='26261fd6-0e6d-4252-a0ea-410b4b4f2eef'
# Abacus_RawData_freezeman-processing_Robot
ABA_EP='fb9e89bd-7bfb-4842-9c4d-9190eb946160'
SRC_BASE_PATH="/lb/robot/research/freezeman-processing"

runfolder=$(jq -r '.run_name' "$run_processing_json")

# Temporary File location, you may want to change it to your scratch for easier clean up.
TRANSFER_JSON_DIR='/lb/project/mugqic/projects/MOH/Transfer_json'
TEMP='/lb/project/mugqic/projects/MOH/TEMP'
TIMESTAMP=$(date +%FT%H.%M.%S)
LOGFILE="${runfolder}_${TIMESTAMP}_${destination}_transfer.log"
LISTFILE="${runfolder}_${TIMESTAMP}_${destination}_transfer.list"
LISTFILE_FULLPATH="${runfolder}_${TIMESTAMP}_${destination}_transfer_fullpath.list"
timestamp_start=$(date "+%Y-%m-%dT%H.%M.%S")

touch "$TEMP/$LOGFILE"
touch "$TEMP/$LISTFILE"
echo "Log file of transfer from Abacus to $destination" > "$TEMP/$LOGFILE"

jq -r '
  .specimen[]?.sample[]? | 
  .sample_name as $sample_name | 
  .readset[]? | 
  .readset_lane as $readset_lane | 
  .file[]? | 
  select(.location_uri != null) | 
  "\($sample_name) \($readset_lane) \(.location_uri | sub("^abacus://"; ""))"
' "$run_processing_json" | while read -r sample_name readset_lane file; do
  file_basename=$(basename "$file")
  file_without_extension="${file_basename%%.sorted.*}"
  file_extension="${file_basename##*.}"
  new_filename="${file_without_extension}_L00${readset_lane}.sorted.${file_extension}"
  if [[ ($file == *.bam || $file == *.bai) && "$file" != *"_L00${readset_lane}.sorted.${file_extension}" ]]; then
    echo "${file#"$SRC_BASE_PATH"} ${DEST_LOC#"$DEST_BASE_PATH"}/$sample_name/$new_filename" >> "$TEMP/$LISTFILE"
    echo "$file $DEST_LOC/$sample_name/$new_filename" >> "$TEMP/$LISTFILE_FULLPATH"
    echo "$file,$sample_name/$new_filename" >> "$TEMP/$LOGFILE"
  else
    echo "${file#"$SRC_BASE_PATH"} ${DEST_LOC#"$DEST_BASE_PATH"}/$sample_name/$file_basename" >> "$TEMP/$LISTFILE"
    echo "$file $DEST_LOC/$sample_name/$file_basename" >> "$TEMP/$LISTFILE_FULLPATH"
    echo "$file,$sample_name/$file_basename" >> "$TEMP/$LOGFILE"
  fi
done

# echo "$TEMP/$LOGFILE $DEST_LOG_LOC/$LOGFILE" >> "$TEMP/$LISTFILE"

if [[ $destination != Abacus ]]; then
  # Load globus module
  module load mugqic/globus-cli/3.24.0
  # Now using Robot account to do the transfer
  ENV_DIR="$HOME/.config/globus_cli"

  case "$destination" in
    Cardinal)
      ENV_FILE="$ENV_DIR/Abacus_to_Cardinal.sh"
      ;;
    Rorqual)
      ENV_FILE="$ENV_DIR/Abacus_to_Rorqual.sh"
      ;;
    Narval)
      ENV_FILE="$ENV_DIR/Abacus_to_Narval.sh"
      ;;
    *)
      echo "ERROR: No Globus env file defined for destination '$destination'." >&2
      exit 1
      ;;
  esac

  if [[ ! -f "$ENV_FILE" ]]; then
    echo "ERROR: Expected environment file '$ENV_FILE' not found." >&2
    exit 1
  fi
  # shellcheck disable=SC1090
  source "$ENV_FILE"

  # Generate and store a UUID for the submission-id
  sub_id="$(globus task generate-submission-id)"
  # Start the batch transfer
  # shellcheck disable=SC2086
  task_id="$(globus transfer --sync-level mtime --jmespath 'task_id' --format=UNIX --submission-id "$sub_id" --label "$runfolder" --batch "$TEMP/$LISTFILE" $ABA_EP $DEST_EP)"

  echo -e "Waiting on 'globus transfer' task '$task_id'.\nTo monitor the transfer see: https://app.globus.org/activity/$task_id/overview"
  globus task wait "$task_id" --polling-interval 60 -H
  # shellcheck disable=SC2181
  if [ $? -eq 0 ]; then
    TRANSFER_JSON="$TRANSFER_JSON_DIR/${LISTFILE/.list/.json}"
    module unload mugqic/globus-cli/3.24.0
    timestamp_end=$(date "+%Y-%m-%dT%H.%M.%S")
    # shellcheck disable=SC2086
    ~/moh_automation/transfer2json.py --input $TEMP/$LISTFILE_FULLPATH --source "abacus" --destination $destination --output $TRANSFER_JSON --operation_cmd_line "globus transfer --sync-level mtime --jmespath 'task_id' --format=UNIX --submission-id $sub_id --label $runfolder --batch $TEMP/$LISTFILE $ABA_EP $DEST_EP" --start $timestamp_start --stop $timestamp_end
    echo "Ingesting transfer $TRANSFER_JSON..."
    # shellcheck disable=SC2086
    pt-cli ingest transfer --input-json $TRANSFER_JSON
  else
    echo "$task_id failed!"
  fi
else
  while IFS= read -r line; do
    # shellcheck disable=SC2046,SC2086
    mkdir -p $(dirname $(echo $line | awk '{print $2}'))
    # shellcheck disable=SC2086
    ln -sf $line
  done < "$TEMP/$LISTFILE"
fi
