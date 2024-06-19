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
  # shellcheck disable=SC2086
  $MOH_path/moh_automation/moh_automation_main/genpipes_deliverables_metrics.py -i $1 -o ${1/.json/_tagged.json}
  module unload mugqic/python/3.11.1
}

genpipes_ingesting() {
  echo "-> Ingesting GenPipes json..."
  # shellcheck disable=SC1091,SC2086
  source $MOH_path/project_tracking_cli/venv/bin/activate
  # shellcheck disable=SC2086
  pt-cli ingest genpipes --input-json $1
  deactivate
}

genpipes_transfer() {
  genpipes_run=$(dirname "$readset_file" | cut -d "/" -f7)
  transfer_log=$(dirname "$readset_file")/transfer.log
  echo "-> Transferring GenPipes run $genpipes_run..."
  {
    # shellcheck disable=SC2086
    (sleep 1 && $MOH_path/moh_automation/moh_automation_main/transfer_GenPipes.sh -r $1 -p $2 -t $3 2>&1) & echo -n "PID: "
    echo $!
    echo "LOG: "
  } >> "$transfer_log"
  echo "-> To follow transfer status see $transfer_log"
}

while getopts 'hc:j:r:l:' OPTION; do
  case "$OPTION" in
  c)
    cluster="$OPTARG"
      if [[ $cluster == abacus ]]; then
        MOH_path="/lb/project/mugqic/projects/MOH"
        if [ -z "${MUGQIC_INSTALL_HOME_DEV:-}" ]; then
          export MUGQIC_INSTALL_HOME_DEV=/lb/project/mugqic/analyste_dev
        fi
        if [ -z "${MUGQIC_INSTALL_HOME_PRIVATE:-}" ]; then
          export MUGQIC_INSTALL_HOME_PRIVATE=/lb/project/mugqic/analyste_private
        fi
      elif [[ $cluster == beluga ]]; then
        MOH_path="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING"
        if [ -z "${MUGQIC_INSTALL_HOME_DEV:-}" ]; then
          export MUGQIC_INSTALL_HOME_DEV=/project/6007512/C3G/analyste_dev
        fi
      elif [[ $cluster == cardinal ]]; then
        MOH_path="/project/def-c3g/MOH"
        if [ -z "${MUGQIC_INSTALL_HOME_DEV:-}" ]; then
          export MUGQIC_INSTALL_HOME_DEV=/project/def-c3g/analyste_dev
        fi
      else
        echo -e "ERROR: Invalid cluster: '$cluster'. It has to be either 'abacus', 'beluga' or 'cardinal'\n"
        usage
      fi
    ;;
  j)
    genpipes_json=$(readlink -f "$OPTARG")
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

if [ -z "${MUGQIC_INSTALL_HOME:-}" ]; then
  export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
fi

if [ -z "${PORTAL_OUTPUT_DIR:-}" ]; then
  export PORTAL_OUTPUT_DIR=$MUGQIC_INSTALL_HOME_DEV/portal_out_dir
fi

module avail 2>&1 | grep -m 1 -q "mugqic"; greprc=$?
if ! [[ $greprc -eq 0 ]]; then
  module use "$MUGQIC_INSTALL_HOME/modulefiles" "$MUGQIC_INSTALL_HOME_DEV/modulefiles"
fi

if [ -z "${JOB_MAIL:-}" ]; then
  export JOB_MAIL=c3g-processing@fakeemail.ca
fi

# Load globus module
module load mugqic/globus-cli/3.24.0
globus_logged=$(globus whoami 2>&1)
if [[ $globus_logged == *"MissingLoginError"* ]]; then
  echo "ERROR: Globus not logged in. Please run 'globus login' in $cluster under robot user. Exiting..."
  exit 1
fi
module unload mugqic/globus-cli/3.24.0

operation_cmd_line=$(jq '.operation_cmd_line' "$genpipes_json")
pipeline=$(echo "$operation_cmd_line" | cut -d' ' -f1 | rev | cut -d'/' -f2 | rev)
protocol=$(echo "$genpipes_json" | cut -d'.' -f2 |  cut -d'_' -f1)

# MAIN folder location
MOH_MAIN="$MOH_path/MAIN"

genpipes_submission_folder=$(dirname "$readset_file")

echo "-> Checking $genpipes_submission_folder..."

module load mugqic/genpipes
if [[ $cluster == beluga ]] || [[ $cluster == cardinal ]] ; then
  log_report_file="${job_list}.tsv"
  # shellcheck disable=SC2046,SC2086
  log_report_output=$(log_report.py $(readlink -f $job_list) --tsv $log_report_file 2>&1)
  failure=$(awk -F'\t' 'NR>1 {print $5"\n"$6"\n"$7}' "$log_report_file" | sort | uniq)
  chmod 660 "$log_report_file"
elif [[ $cluster == abacus ]]; then
  log_report_file="${job_list}.txt"
  # shellcheck disable=SC2086
  log_report_output=$(log_report.pl -nos $job_list)
  failure=$(echo "$log_report_output" | grep -v "^#" | awk -F'\t' '{print $5}' | sort | uniq)
  echo "$log_report_output" > "$MOH_MAIN/job_output/${job_list}.txt"
  chmod 660 "$MOH_MAIN/job_output/${job_list}.txt"
fi
# echo "failure: $failure"
if [[ $failure == *"FAILED"* ]] || [[ $failure == *"TIMEOUT"* ]]; then
  echo "WARNING: Failure found in $job_list Cf. $log_report_file"
  # Let's tag GenPipes + Ingest GenPipes
  genpipes_tagging "$genpipes_json"
  genpipes_ingesting "${genpipes_json/.json/_tagged.json}"
  touch "${genpipes_submission_folder}.checked"
  chmod 660 "${genpipes_submission_folder}.checked"
elif [[ $failure == *"ACTIVE"* ]] || [[ $failure == *"RUNNING"* ]] || [[ $failure == *"PENDING"* ]]; then
  # Let's skip and wait
  echo "INFO: Job(s) still running Cf. $log_report_file"
elif [[ -z $failure ]] || [[ $failure == *"COMPLETED"* ]]; then
  # Let's tag GenPipes + Ingest GenPipes
  genpipes_tagging "$genpipes_json"
  genpipes_ingesting "${genpipes_json/.json/_tagged.json}"
  # Let's transfer GenPipes only if NOT on beluga
  if ! [[ $cluster == beluga ]]; then
    genpipes_transfer "$readset_file" "$pipeline" "$protocol"
  fi
  touch "${genpipes_submission_folder}.checked"
  chmod 660 "${genpipes_submission_folder}.checked"
else
  echo "ERROR: Unknown status Cf. $log_report_file"
fi
