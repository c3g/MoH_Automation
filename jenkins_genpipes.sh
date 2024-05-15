#!/bin/bash
set -eu -o pipefail


usage() {
  echo "script usage: jenkins_genpipes.sh -h [-c cluster] [-p pipeline] [-t protocol] [-i input_file]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -c <cluster>                     Cluster name on which the analysis will be run."
  echo " -p <pipeline>                    Pipeline name to be used for the analysis."
  echo " -t <protocol>                    Protocol to be used for the analysis. (Optional)"
  echo " -i <input_file>                  Path to Input File to be used for the analysis. This file is a csv file with 1st
                                          column being the readset file and the 2nd (optional) column being the pair file."
  exit 1
  }

while getopts 'hc:p::t:i:' OPTION; do
  case "$OPTION" in
    c)
      cluster="$OPTARG"
      if [[ $cluster == beluga ]]; then
        path="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN"
        # path="/scratch/stretenp/tmp/MoH_GenPipes"
        scheduler="slurm"
        if [ -z "$MUGQIC_INSTALL_HOME_DEV" ]; then
          export MUGQIC_INSTALL_HOME_DEV=/project/6007512/C3G/analyste_dev
        fi

      elif [[ $cluster == abacus ]]; then
        # path="/lb/project/mugqic/projects/MOH/MAIN"
        path="/lb/scratch/pstretenowich/MOH/MAIN"
        scheduler="pbs"
        if [ -z "$MUGQIC_INSTALL_HOME_DEV" ]; then
          export MUGQIC_INSTALL_HOME_DEV=/lb/project/mugqic/analyste_dev
        fi
        if [ -z "$MUGQIC_INSTALL_HOME_PRIVATE" ]; then
          export MUGQIC_INSTALL_HOME_PRIVATE=/lb/project/mugqic/analyste_private
        fi
      else
        echo -e "ERROR: Invalid cluster: '$cluster'.\n"
        usage
      fi
      # echo "cluster: $cluster"
      ;;
    p)
      pipeline="$OPTARG"
      case "$pipeline" in
        tumor_pair | rnaseq | rnaseq_light)
          # echo "pipeline: $pipeline"
          ;;
        *)
          echo -e "ERROR: Invalid pipeline: '$pipeline'.\n"
          usage
          ;;
      esac
      ;;
    t)
      protocol="$OPTARG"
      case "$protocol" in
        ensemble | sv | cancer)
          # echo "protocol: $protocol"
          ;;
        *)
          echo -e "ERROR: Invalid protocol: '$protocol'.\n"
          usage
          ;;
      esac
      ;;
    i)
      input_file=$(realpath "$OPTARG")
      # echo "input_file: $input_file"
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
if [ -z "$cluster" ] || [ -z "$pipeline" ] || [ -z "$input_file" ]; then
  echo -e "ERROR: Missing mandatory arguments -c and -p and -i.\n"
  usage
fi

# wrong protocol with wrong pipeline
if [ "$pipeline" = "tumor_pair" ] && ! [[ "$protocol" == "ensemble" || "$protocol" == "sv" ]]; then
  echo -e "ERROR: pipeline: '$pipeline' only accepts protocol: 'ensemble' or 'sv', '$protocol' provided.\n"
  usage
elif [ "$pipeline" = "rnaseq" ] && [ "$protocol" != "cancer" ]; then
  echo -e "ERROR: pipeline: '$pipeline' only accepts protocol: 'cancer', '$protocol' provided.\n"
  usage
fi

if [ -z "$MUGQIC_INSTALL_HOME" ]; then
  export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
fi

if [ -z "$PORTAL_OUTPUT_DIR" ]; then
  export PORTAL_OUTPUT_DIR=$MUGQIC_INSTALL_HOME_DEV/portal_out_dir
fi

module avail 2>&1 | grep -m 1 -q "mugqic"; greprc=$?
if ! [[ $greprc -eq 0 ]]; then
  module use "$MUGQIC_INSTALL_HOME/modulefiles" "$MUGQIC_INSTALL_HOME_DEV/modulefiles"
fi

if [ -z "$JOB_MAIL" ]; then
  export JOB_MAIL=c3g-processing@fakeemail.ca
fi

##################################################
# Initialization
chunk_submit=false
module purge
module load mugqic/python/3.10.2


export MUGQIC_PIPELINES_HOME=${path}/genpipes_moh/genpipes

if [[ $cluster == beluga ]]; then
  beluga_ini="$MUGQIC_PIPELINES_HOME/pipelines/common_ini/beluga.ini"
elif [[ $cluster == abacus ]]; then
  beluga_ini=""
fi

cd "$path"

while IFS=, read -r readset_file pair_file; do
  timestamp=$(date "+%Y-%m-%dT%H.%M.%S")
  timestamp_find_format=$(date "+%Y-%m-%d %H:%M:%S")
  patient=$(awk 'NR==2, match($1, /^((MoHQ-(JG|HM|CM|GC|MU|MR|IQ|XX)-\w+)-\w+)/) {print substr($1, RSTART, RLENGTH)}' "$readset_file")
  # sample=$(awk 'NR>1{print $1}' "$readset_file")
  # echo $sample
  echo "-> Running GenPipes for ${patient}..."
  # GenPipes call
  if test "$pipeline" == rnaseq_light; then
    pipeline_name=RnaSeqLight
    link_folder="${path}/genpipes_submission/${pipeline_name}.${timestamp}"
    mkdir -p "$link_folder"
    # rnaseq_light
    genpipes_file=RnaSeqLight_${patient}.${timestamp}.sh
    # shellcheck disable=SC2086
    $MUGQIC_PIPELINES_HOME/pipelines/rnaseq_light/rnaseq_light.py \
-s 1-4 \
-c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq_light/rnaseq_light.base.ini $beluga_ini \
$MUGQIC_PIPELINES_HOME/pipelines/common_ini/Homo_sapiens.GRCh38.ini \
RNA_light.custom.ini \
-j $scheduler \
-r $readset_file \
-g $genpipes_file \
--json-pt &> "${path}/genpipes_logs/${patient}.${timestamp}.log"
    chunk_submit=true
    after_genpipes_call_timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    json_regex="${path}/json/${pipeline_name}_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].json"
    job_list_regex="${path}/job_output/${pipeline_name}.job_list.[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9]"
    trace_ini_regex="${path}/${pipeline_name}.[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].config.trace.ini"
  elif test "$pipeline" == rnaseq; then
    pipeline_name=RnaSeq
    link_folder="${path}/genpipes_submission/${pipeline_name}.${protocol}.${timestamp}"
    mkdir -p "$link_folder"
    # rnaseq
    genpipes_file=RnaSeq.${protocol}_${patient}.${timestamp}.sh
    # shellcheck disable=SC2086
    $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.py \
-t $protocol \
-s 1-8,11-28 \
-c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini $beluga_ini \
$MUGQIC_PIPELINES_HOME/pipelines/common_ini/Homo_sapiens.GRCh38.ini \
RNA_cancer.custom.ini \
-j $scheduler \
-r $readset_file \
-g $genpipes_file \
--json-pt &> "${path}/genpipes_logs/${patient}.${timestamp}.log"
    chunk_submit=true
    after_genpipes_call_timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    json_regex="${path}/json/${pipeline_name}.${protocol}_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].json"
    job_list_regex="${path}/job_output/${pipeline_name}.${protocol}.job_list.[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9]"
    trace_ini_regex="${path}/${pipeline_name}.${protocol}.[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].config.trace.ini"
  elif test "$pipeline" == tumor_pair; then
    pipeline_name=TumorPair
    link_folder="${path}/genpipes_submission/${pipeline_name}.${protocol}.${timestamp}"
    mkdir -p "$link_folder"
    # tumor_pair
    if test "$protocol" == ensemble; then
      steps="5-13,15-38"
      custom_ini="TP_ensemble.custom.ini"
    elif test "$protocol" == sv; then
      steps="12-16"
      custom_ini="TP_sv.custom.ini"
    fi
    genpipes_file=TumorPair.${protocol}_${patient}.${timestamp}.sh
    # shellcheck disable=SC2086
    $MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.py \
-t $protocol \
-s $steps \
-c $MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.base.ini \
$MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.extras.ini $beluga_ini \
$MUGQIC_PIPELINES_HOME/pipelines/common_ini/Homo_sapiens.GRCh38.ini \
$custom_ini \
-j $scheduler \
-r $readset_file \
-p $pair_file \
-g $genpipes_file \
--json-pt &> "${path}/genpipes_logs/${patient}.${timestamp}.log"
    chunk_submit=true
    after_genpipes_call_timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    json_regex="${path}/json/${pipeline_name}.${protocol}_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].json"
    job_list_regex="${path}/job_output/${pipeline_name}.${protocol}.job_list.[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9]"
    trace_ini_regex="${path}/${pipeline_name}.${protocol}.[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].config.trace.ini"
  fi
  # Chunking & Submission
  if test $chunk_submit == true && test -f "$genpipes_file"; then
    today=$(date "+%Y-%m-%dT")
    chmod 664 -- *."$protocol"."$today"*.config.trace.ini
    chmod 774 "$genpipes_file"
    echo "-> Chunking GenPipes for ${patient}..."
    "$MUGQIC_PIPELINES_HOME"/utils/chunk_genpipes.sh "$genpipes_file" "${patient}.${timestamp}_chunks" &> "${path}/genpipes_logs/${patient}.${timestamp}_chunks.log"
    chmod 775 "${patient}.${timestamp}_chunks"
    chmod 664 "${patient}.${timestamp}_chunks"/*
    echo "-> Submitting GenPipes for ${patient}..."
    cat /dev/null > "${path}/genpipes_logs/${patient}.${timestamp}.txt"
    {
      (sleep 1 && "$MUGQIC_PIPELINES_HOME"/utils/submit_genpipes "${patient}.${timestamp}_chunks" 2>&1) & echo -n "PID: "
      echo $!
      echo "PATIENT: ${patient}"
      echo "LOG: "
    } >> "${path}/genpipes_logs/${patient}.${timestamp}.txt"
    chmod 664 "${patient}.${timestamp}.txt"
    # Finding GenPipes json and symklink into genpipes_submission along with the readset file
    maybe_json=$(find "${path}/json" -type f -regex "$json_regex" -newermt "$timestamp_find_format" | sort | tail -n 1)
    ln -s "$readset_file" "$link_folder"/.
    ln -s "$maybe_json" "$link_folder"/.
    # Need to wait for the scheduler to return the job IDs and so have the job_list file generated
    echo "-> Waiting for job_list to be created for ${patient}..."
    maybe_job_list=""
    while [ -z "$maybe_job_list" ]; do
      sleep 1
      maybe_job_list=$(find "${path}/job_output" -type f -regex "$job_list_regex" -newermt "$after_genpipes_call_timestamp" | sort | tail -n 1)
    done
    ln -s "$maybe_job_list" "$link_folder"/.
    # Do some cleaning
    mv "$genpipes_file" "${path}/genpipes_files"
    maybe_trace_ini=$(find "${path}" -type f -regex "$trace_ini_regex" -newermt "$timestamp_find_format" | sort | tail -n 1)
    mv "$maybe_trace_ini" "${path}/genpipes_inis"
    mv "${patient}.${timestamp}_chunks" "${path}/genpipes_logs"
  fi
done < "${input_file}"
