#!/bin/bash
set -eu -o pipefail


usage() {
  echo "script usage: ProjectTracking_RunProcessing.sh -h [-r runfolder] [-l lane] [-s sample] [-x xsample]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -r <runfolder>                   RunFolder found in abacus under /lb/robot/research/freezeman-processing/<sequencer>/<year>/<runfolder>."
  echo " -d <destination>                 Destination of the transfer (Abacus, Beluga or Cardinal). For now it will also always be transferred to Beluga for Old DB compatibility."
  echo " -l <lane>                        Lane(s) to be ingested (default: all)."
  echo " -s <sample>                      Sample Name(s) (as they appear in the json file from Freezeman) (default: all)."
  echo " -x <xsample>                     Sample Name(s) to be EXCLUDED (as they appear in the json file from Freezeman) (default: none)."
  exit 1
  }

while getopts 'hr:d::l::s::x:' OPTION; do
  case "$OPTION" in
    r)
      runfolder="$OPTARG"
      ;;
    d)
      destination="$OPTARG"
      ;;
    l)
      lane+=("$OPTARG")
      ;;
    s)
      sample+=("$OPTARG")
      ;;
    x)
      xsample+=("$OPTARG")
      ;;
    h)
      usage
      ;;
    ?)
      usage
      ;;
  esac
done

# location on Beluga. CURRENTLY VERY IMPORTANT. DO NOT CHANGE OR IT WILL BREAK THE DATABASE
# SERIOUSLY DON'T CHANGE IT.
# PLEASE DONT.
if [[ $destination = Beluga ]]; then
    # Beluga main folder location
    DEST_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/raw_reads"
    # Beluga log file location
    DEST_LOG_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer"
    # Beluga Endpoint
    DEST_EP='278b9bfe-24da-11e9-9fa2-0a06afd4a22e'
elif [[ $destination = Cardinal ]]; then
    # Cardinal main folder location
    DEST_LOC="/project/60007/MOH/MAIN/raw_reads"
    # Cardinal log file location
    DEST_LOG_LOC="/project/60007/MOH/log_files/transfer"
    # Cardinal Endpoint
    DEST_EP='26f926d9-6216-4e84-9037-a5c9567b5707'
elif [[ $destination = Abacus ]]; then
    # Abacus main folder location
    DEST_LOC="/lb/project/mugqic/projects/MOH/MAIN/raw_reads"
    # Abacus log file location
    DEST_LOG_LOC="/lb/project/mugqic/projects/MOH/log_files/transfer"
    # Abacus Endpoint
    DEST_EP=''
fi

run_processing2json_args=""
if declare -p lane >&/dev/null; then
  run_processing2json_args="${run_processing2json_args} -l ${lane[*]}"
fi
if declare -p sample >&/dev/null; then
  run_processing2json_args="${run_processing2json_args} -s ${sample[*]}"
fi
if declare -p xsample >&/dev/null; then
  run_processing2json_args="${run_processing2json_args} -x ${xsample[*]}"
fi

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
export MUGQIC_INSTALL_HOME_DEV=/lb/project/mugqic/analyste_dev
module use "$MUGQIC_INSTALL_HOME/modulefiles" "$MUGQIC_INSTALL_HOME_DEV/modulefiles"

##################################################
# Initialization
module purge
# module load mugqic/python/3.10.2
# Using python from pt_cli env
# shellcheck disable=SC1090
source ~/project_tracking_cli/venv/bin/activate

path=/lb/project/mugqic/projects/MOH/RunProcessing

cd "$path"

## Prepare run list
folder_prefix=/lb/robot/research/freezeman-processing
input=$(find "$folder_prefix"/*/*/ -type d -name "$runfolder")
echo "-> Processing $runfolder..."
# Temporary File location, you may want to change it to your scratch for easier clean up.
TEMP='/lb/project/mugqic/projects/MOH/TEMP'
TIMESTAMP=$(date +%FT%H.%M.%S)
LOGFILE="${runfolder}_${TIMESTAMP}_transfer.log"
LISTFILE="${runfolder}_${TIMESTAMP}_transfer.list"

# abacus Endpoint
ABA_EP='26261fd6-0e6d-4252-a0ea-410b4b4f2eef'

if [ -s "$input" ]; then
  run_processing_json=$path/$runfolder.json
  # Json creation from run csv file
  # shellcheck disable=SC2086
  ~/moh_automation/fms_run_processing2json.py $run_processing2json_args --input $input --output $run_processing_json
  chmod 664 "$run_processing_json"
  # Using client to add new runs to database
  # shellcheck disable=SC2086
  ret="$(pt-cli ingest run_processing --input-json $run_processing_json 2>&1 || true)"
  echo -e "$ret"
  if ! [[ $ret == *"has to be unique"* ]] && [[ $ret == *"BadRequestError"* ]]; then
    exit 1
  fi
else
  echo "--> ERROR: Missing folder $folder_prefix/*/*/$runfolder"
  exit 1
fi

touch "$TEMP/$LOGFILE"
touch "$TEMP/$LISTFILE"
echo "Log file of transfer from Abacus to $destination" > "$TEMP/$LOGFILE"
echo "Transferred From $input" >> "$TEMP/$LOGFILE"

for sample_name in $(jq -r '.. | .sample_name? // empty' $run_processing_json); do
  for file in $(jq -r --arg sample_name "$sample_name" '.. | select(.sample_name? == $sample_name) | .readset[]?.file[]?.location_uri? // empty | sub("^abacus://"; "")' $run_processing_json); do
    echo "$file $DEST_LOC/$sample_name/$(basename $file)" >> "$TEMP/$LISTFILE"
    echo "$file,$sample_name/$(basename $file)" >> "$TEMP/$LOGFILE"
  done
done

echo "$TEMP/$LOGFILE $DEST_LOG_LOC/$LOGFILE" >> "$TEMP/$LISTFILE"

if [[ $destination != Abacus ]]; then
  # Load globus module
  module load mugqic/globus-cli/3.24.0

  # Generate and store a UUID for the submission-id
  sub_id="$(globus task generate-submission-id)"

  # Start the batch transfer
  task_id="$(globus transfer --sync-level mtime --jmespath 'task_id' --format=UNIX --submission-id "$sub_id" --label "$runfolder" --batch "$TEMP/$LISTFILE" $ABA_EP $DEST_EP)"

  echo "Waiting on 'globus transfer' task '$task_id'"
  globus task wait "$task_id" --polling-interval 60 -H
  # shellcheck disable=SC2181
  if [ $? -eq 0 ]; then
    module unload mugqic/globus-cli/3.24.0
    # shellcheck disable=SC2086
    ~/moh_automation/transfer2json.py --input $TEMP/$LISTFILE --source "abacus" --destination $destination --output /lb/project/mugqic/projects/MOH/Transfer_json/${LISTFILE/.list/.json} --operation_cmd_line "globus transfer --sync-level mtime --jmespath 'task_id' --format=UNIX --submission-id $sub_id --label $runfolder --batch $TEMP/$LISTFILE $ABA_EP $DEST_EP"
    # shellcheck disable=SC2086
    pt-cli ingest transfer --input-json /lb/project/mugqic/projects/MOH/Transfer_json/${LISTFILE/.list/.json}
  else
    echo "$task_id failed!"
  fi
else
  while IFS= read -r line; do
    mkdir -p $(dirname $(echo $line | awk '{print $2}'))
    ln -sf $line
  done < "$TEMP/$LISTFILE"
fi
