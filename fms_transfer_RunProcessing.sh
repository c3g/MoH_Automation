#!/usr/bin/env bash

THIS_SCRIPT=$(basename "$0")

usage() {
  echo "script usage: $THIS_SCRIPT -h [-l run_processing_json] [-d destination] [-s sample]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -r <run_processing_json>         Run Processing json."
  echo " -d <destination>                 Destination for the transfer (either Beluga or Cardinal or Abacus)."
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

if ! [[ $destination =~ Cardinal|Beluga|Abacus ]]; then
    echo -e "ERROR: Invalid destination: '$destination'. It has to be either Beluga, Cardinal or Abacus.\n"
    usage
fi

# location on Beluga. CURRENTLY VERY IMPORTANT. DO NOT CHANGE OR IT WILL BREAK THE DATABASE
# SERIOUSLY DON'T CHANGE IT.
# PLEASE DONT.
if [[ $destination = Beluga ]]; then
    # Beluga Endpoint
    DEST_EP='278b9bfe-24da-11e9-9fa2-0a06afd4a22e'
    # Beluga_MoH_Robot Endpoint
    # DEST_EP='43c48ed2-8c4e-4c4d-bd4d-f29ace1c1a5e'
    # Beluga base path
    DEST_BASE_PATH="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING"
    # Beluga main folder location
    DEST_LOC="${DEST_BASE_PATH}/MAIN/raw_reads"
    # Beluga log file location
    DEST_LOG_LOC="${DEST_BASE_PATH}/DATABASE/log_files/transfer"
elif [[ $destination = Cardinal ]]; then
    # Cardinal Endpoint
    DEST_EP='a6df16bc-4e7f-4784-9afa-8ceb7b20b7c0'
    # Cardinal_MoH_Robot Endpoint
    # DEST_EP='995863d1-97fe-4af2-bc48-c9bf91c33f08'
    # Cardinal base path
    DEST_BASE_PATH="/project/60007/MOH"
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

# abacus Endpoint
ABA_EP='26261fd6-0e6d-4252-a0ea-410b4b4f2eef'

runfolder=$(jq -r '.run_name' "$run_processing_json")

# Temporary File location, you may want to change it to your scratch for easier clean up.
TEMP='/lb/project/mugqic/projects/MOH/TEMP'
TIMESTAMP=$(date +%FT%H.%M.%S)
LOGFILE="${runfolder}_${TIMESTAMP}_${destination}_transfer.log"
LISTFILE="${runfolder}_${TIMESTAMP}_${destination}_transfer.list"

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
    echo "$file $DEST_LOC/$sample_name/$new_filename" >> "$TEMP/$LISTFILE"
    echo "$file,$sample_name/$new_filename" >> "$TEMP/$LOGFILE"
  else
    echo "$file $DEST_LOC/$sample_name/$file_basename" >> "$TEMP/$LISTFILE"
    echo "$file,$sample_name/$file_basename" >> "$TEMP/$LOGFILE"
  fi
done

echo "$TEMP/$LOGFILE $DEST_LOG_LOC/$LOGFILE" >> "$TEMP/$LISTFILE"

if [[ $destination != Abacus ]]; then
  # Load globus module
  module load mugqic/globus-cli/3.24.0

  # Generate and store a UUID for the submission-id
  sub_id="$(globus task generate-submission-id)"
  # Start the batch transfer
  # shellcheck disable=SC2086
  task_id="$(globus transfer --sync-level mtime --jmespath 'task_id' --format=UNIX --submission-id "$sub_id" --label "$runfolder" --batch "$TEMP/$LISTFILE" $ABA_EP $DEST_EP)"

  echo -e "Waiting on 'globus transfer' task '$task_id'.\nTo monitor the transfer see: https://app.globus.org/activity/$task_id/overview"
  globus task wait "$task_id" --polling-interval 60 -H
  # shellcheck disable=SC2181
  if [ $? -eq 0 ]; then
    module unload mugqic/globus-cli/3.24.0
    # shellcheck disable=SC2086
    ~/moh_automation/transfer2json.py --input $TEMP/$LISTFILE --source "abacus" --destination $destination --output /lb/project/mugqic/projects/MOH/Transfer_json/${LISTFILE/.list/.json} --operation_cmd_line "globus transfer --sync-level mtime --jmespath 'task_id' --format=UNIX --submission-id $sub_id --label $runfolder --batch $TEMP/$LISTFILE $ABA_EP $DEST_EP"
    echo "Ingesting transfer /lb/project/mugqic/projects/MOH/Transfer_json/${LISTFILE/.list/.json}..."
    # shellcheck disable=SC2086
    pt-cli ingest transfer --input-json /lb/project/mugqic/projects/MOH/Transfer_json/${LISTFILE/.list/.json}
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
