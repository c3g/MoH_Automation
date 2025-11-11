#!/bin/bash
set -e -o pipefail

usage() {
  echo "script usage: ProjectTracking_RunProcessing.sh -h [-r runfolder] [-l lane] [-s sample] [-x xsample]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -r <runfolder>                   RunFolder found in abacus under /lb/robot/research/freezeman-processing/<sequencer>/<year>/<runfolder>."
  echo " -d <destination>                 Destination of the transfer (Abacus, Cardinal or Rorqual)."
  echo " -n <nucleic_acid_type>           nucleic_acid_type to be considered for the transfer (either DNA or RNA, Default: ALL)."
  echo " -l <lane>                        Lane(s) to be ingested (default: all)."
  echo " -s <sample>                      Sample Name(s) (as they appear in the json file from Freezeman) (default: all)."
  echo " -x <xsample>                     Sample Name(s) to be EXCLUDED (as they appear in the json file from Freezeman) (default: none)."
  exit 1
  }

while getopts 'hr:d::n::l::s::x:' OPTION; do
  case "$OPTION" in
    r)
      runfolder="$OPTARG"
      ;;
    d)
      destination="$OPTARG"
      ;;
    n)
      nucleic_acid_type="$OPTARG"
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
input=$(find "$folder_prefix"/*/*/ -maxdepth 1 -type d -name "$runfolder")
echo "-> Processing $runfolder..."

TIMESTAMP=$(date +%FT%H.%M.%S)

if [ -s "$input" ]; then
  run_processing_json=$path/$runfolder.json
  # Json creation from run csv file
  if [ -z "$nucleic_acid_type" ]; then
    # shellcheck disable=SC2086
    ~/moh_automation/fms_run_processing2json.py $run_processing2json_args --input $input --output $run_processing_json
  else
    # shellcheck disable=SC2086
    ~/moh_automation/fms_run_processing2json.py $run_processing2json_args --input $input --output $run_processing_json --nucleic_acid_type "$nucleic_acid_type"
  fi
  
  chmod 664 "$run_processing_json"
  # Sample check before ingestion: make sure sample name is correct aka an existing one is not renamed
  ~/moh_automation/pt_check_sample.sh -i $run_processing_json
  # Using client to add new runs to database
  # shellcheck disable=SC2086
  ret="$(pt-cli ingest run_processing --input-json $run_processing_json 2>&1 || true)"
  echo -e "$ret" | tee "${run_processing_json/.json/_${TIMESTAMP}_run_processing_ingestion.log}"
  if ! [[ $ret == *"has to be unique"* ]] && [[ $ret == *"BadRequestError"* ]]; then
    exit 1
  fi
else
  echo "--> ERROR: Missing folder $folder_prefix/*/*/$runfolder"
  exit 1
fi

nohup ~/moh_automation/fms_transfer_RunProcessing.sh -r "$run_processing_json" -d "$destination" > "${run_processing_json/.json/_${destination}_${TIMESTAMP}_transfer.log}" 2>&1 &
echo "Transfer started towards $destination. See log file ${run_processing_json/.json/_${destination}_${TIMESTAMP}_transfer.log}"