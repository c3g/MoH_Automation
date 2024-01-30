#!/bin/bash
set -eu -o pipefail


usage() {
  echo "script usage: run_processing2json_wrapper.sh -h [-c cluster] [-p pipeline] [-t protocol] [-i input_file]"
  echo "Usage:"
  echo " -h                               Display this help message."
  exit 1
  }

while getopts 'h' OPTION; do
  case "$OPTION" in
    h)
      usage
      ;;
    ?)
      usage
      ;;
  esac
done

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
echo "-> Detecting new runs..."
runs_folder=/lb/robot/research/processing/novaseq
cat /dev/null > all.runs.txt.tmp

find "$runs_folder"/*/ -maxdepth 1 -iname "*M[O\|o]H*Run*" -type d -exec basename {} \; > all.runs.txt.tmp
# Diff between ingested runs and all runs to find new runs not ingested
diff -Bw <(sort ingested.runs.txt) <(sort all.runs.txt.tmp) | grep "^>" | sed s:"> ":: > new.runs.tmp && ret=$? || ret=$?

if ((ret >= 2)); then
 exit "$ret"
fi

if [ -s new.runs.tmp ]; then
  while IFS= read -r run; do
    echo "-> Processing $run"
    input=$(find "$runs_folder"/*/"$run/" -name "$run-run.align_bwa_mem.csv")
    if [ -s "$input" ]; then
      # Json creation from run csv file
      # shellcheck disable=SC2086
      ~/moh_automation/run_processing2json.py --input $input --output $path/$run.json
      chmod 664 "$path/$run.json"
      # Using client to add new runs to database
      # shellcheck disable=SC2086
      pt_cli route /project
      # pt_cli ingest run_processing --input-json $path/$run.json
      # echo "$run" >> ingested.runs.txt
    else
      echo "--> WARNING: Missing $runs_folder/*/$run/$run-run.align_bwa_mem.csv file, skipping..."
    fi
  done < new.runs.tmp
else
  echo "No new runs detected."
fi

rm all.runs.txt.tmp new.runs.tmp
