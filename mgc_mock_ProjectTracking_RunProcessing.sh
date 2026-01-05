#!/bin/bash
set -eu -o pipefail

usage() {
  echo "script usage: ProjectTracking_RunProcessing.sh -h [-r runid] [-l lane] [-s sample] [-x xsample]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -r <runid>                       RunID found in abacus under /lb/robot/research/processing/novaseq."
  echo " -l <lane>                        Lane(s) to be ingested (default: all)."
  echo " -s <sample>                      Sample Name(s) (as they appear in the file <runid>-run.align_bwa_mem.csv) (default: all)."
  echo " -x <xsample>                     Sample Name(s) to be EXCLUDED (as they appear in the file <runid>-run.align_bwa_mem.csv) (default: none)."
  exit 1
  }

while getopts 'hr::l::s::x:' OPTION; do
  case "$OPTION" in
    r)
      runid="$OPTARG"
      echo "runid: $runid"
      ;;
    l)
      lane+=("$OPTARG")
      echo "lane: ${lane[*]}"
      ;;
    s)
      sample+=("$OPTARG")
      echo "sample: ${sample[*]}"
      ;;
    x)
      xsample+=("$OPTARG")
      echo "xsample: ${xsample[*]}"
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
runs_folder=/lb/project/mugqic/projects/MOH/GQ_STAGING/mgc_mock_runs
input=$(find "$runs_folder"/"$runid/" -name "$runid-run.align_bwa_mem.csv")
echo "-> Processing $runid, file $input..."

TIMESTAMP=$(date +%FT%H.%M.%S)

if [ -s "$input" ]; then
  # Json creation from run csv file
  # shellcheck disable=SC2086
  ~/moh_automation/mgc_mock_run_processing2json.py $run_processing2json_args --input $input --output $path/$runid.json
  chmod 664 "$path/$runid.json"
  ~/moh_automation/pt_check_sample.sh $path/$runid.json
  # Using client to add new runs to database
  # shellcheck disable=SC2086
  ret="$(pt-cli ingest run_processing --input-json $path/$runid.json 2>&1 || true)"
  echo -e "$ret" | tee "${path}/${runid}_${TIMESTAMP}_run_processing_ingestion.log"
  if ! [[ $ret == *"has to be unique"* ]] && [[ $ret == *"BadRequestError"* ]]; then
    exit 1
  fi
else
  echo "--> ERROR: Missing $runs_folder/*/$runid/$runid-run.align_bwa_mem.csv file"
  exit 1
fi

nohup ~/moh_automation/mgc_mock_transfer_RunProcessing.sh -r "$run_processing_json" -d "$destination" > "${run_processing_json/.json/_${destination}_${TIMESTAMP}_transfer.log}" 2>&1 &
echo "Transfer started towards $destination. See log file ${run_processing_json/.json/_${destination}_${TIMESTAMP}_transfer.log}"