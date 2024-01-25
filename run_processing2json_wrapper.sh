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
module load mugqic/python/3.10.2

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
    input=$(find "$runs_folder"/*/"$run/" -name "$run-novaseq-run.align_bwa_mem.csv")
    echo "./run_processing2json.py --input $input --output $path/$run.json"
  done < new.runs.tmp
else
  echo "No new runs detected."
fi

rm all.runs.txt.tmp new.runs.tmp