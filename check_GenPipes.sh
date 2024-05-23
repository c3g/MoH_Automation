#!/usr/bin/env bash

THIS_SCRIPT=$(basename "$0")

usage() {
  echo "script usage: $THIS_SCRIPT -h [-c cluster] [-j genpipes_json] [-r readset_file] [-l job_list]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -c <cluster>                     Cluster name on which the check is expected /!\ MANDATORY /!\."
  echo " -j <genpipes_json>               json file generated by GenPipes /!\ MANDATORY /!\."
  echo " -r <readset_file>                Readset File used to submit GenPipes /!\ MANDATORY /!\."
  echo " -l <job_list>                    job_list file generated by GenPipes /!\ MANDATORY /!\."
  exit 1
  }

genpipes_tagging() {
  echo "-> Tagging GenPipes json..."
  module load mugqic/python/3.11.1
  echo "$MOH_path/moh_automation/moh_automation_main/genpipes_deliverables_metrics.py -i $1 -o ${1/.json/_tagged.json}"
  module unload mugqic/python/3.11.1
}

genpipes_ingesting() {
  echo "-> Ingesting GenPipes json..."
  # shellcheck disable=SC1091,SC2086
  source $MOH_path/project_tracking_cli/venv/bin/activate
  echo "pt-cli ingest genpipes --input-json $1"
  deactivate
}

genpipes_transfer() {
  echo "-> Transfering GenPipes..."
  echo "$MOH_path/moh_automation/moh_automation_main/transfer_GenPipes.sh -r $1 -p $2 -t $3"
}

while getopts 'hc:j:r:l:' OPTION; do
  case "$OPTION" in
  c)
    cluster="$OPTARG"
    if [[ $cluster == beluga ]]; then
        MOH_path="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING"
      elif [[ $cluster == abacus ]]; then
        MOH_path="/lb/project/mugqic/projects/MOH"
      elif [[ $cluster == cardinal ]]; then
        MOH_path="/project/def-c3g/MOH"
      else
        echo -e "ERROR: Invalid cluster: '$cluster'.\n"
        usage
      fi
    ;;
  j)
    genpipes_json="$OPTARG"
    ;;
  r)
    readset_file="$OPTARG"
    ;;
  l)
    job_list="$OPTARG"
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
if [ ! "$cluster" ] || [ ! "$readset_file" ] || [ ! "$genpipes_json" ] || [ ! "$job_list" ]; then
  echo -e "ERROR: Missing mandatory arguments -c and/or -r and/or -j and/or -l.\n"
  usage
fi

operation_cmd_line=$(jq '.operation_cmd_line' "$genpipes_json")
pipeline=$(echo "$operation_cmd_line" | cut -d' ' -f1 | rev | cut -d'/' -f2 | rev)
protocol=$(echo "$genpipes_json" | cut -d'.' -f2 |  cut -d'_' -f1)

# MAIN folder location
MOH_MAIN="$MOH_path/MAIN"

module load mugqic/genpipes
if [[ $cluster == beluga ]] || [[ $cluster == cardinal ]] ; then
  log_report_file="${job_list}.tsv"
  # shellcheck disable=SC2086
  log_report_output=$(log_report.py $job_list --tsv $log_report_file)
  failure=$log_report_output
elif [[ $cluster == abacus ]]; then
  log_report_file="${job_list}.txt"
  # shellcheck disable=SC2086
  log_report_output=$(log_report.pl -nos $job_list)
  failure=$(echo "$log_report_output" | grep -v "^#")
  echo "$log_report_output" > "$MOH_MAIN/job_output/${job_list}.txt"
fi
# echo "failure: $failure"
if [[ -z $failure ]]; then
  # Let's tag GenPipes + Ingest GenPipes
  genpipes_tagging "$genpipes_json"
  genpipes_ingesting "${genpipes_json/.json/_tagged.json}"
  # Let's transfer GenPipes only if NOT on beluga
  if ! [[ $cluster == beluga ]]; then
    genpipes_transfer "$readset_file" "$pipeline" "$protocol"
  fi
elif [[ $failure == *"FAILED"* ]] || [[ $failure == *"TIMEOUT"* ]]; then
  # Let's tag GenPipes + Ingest GenPipes
  genpipes_tagging "$genpipes_json"
  genpipes_ingesting "${genpipes_json/.json/_tagged.json}"
  echo "WARNING: Failure found in $job_list Cf. $MOH_MAIN/job_output/$log_report_file"
elif [[ $failure == *"ACTIVE"* ]] || [[ $failure == *"RUNNING"* ]]; then
  # Let's skip and wait
  echo "WARNING: Job(s) still running for $job_list"
else
  echo "ERROR: Unknown status in $job_list"
fi
