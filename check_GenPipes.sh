#!/usr/bin/env bash

THIS_SCRIPT=$(basename "$0")

usage() {
  echo "script usage: $THIS_SCRIPT -h [-j genpipes_json] [-r readset_file]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -r <readset_file>                Readset File used to submit GenPipes /!\ MANDATORY /!\."
  echo " -p <pipeline>                    GenPipes pipeline to be transfered (either: tumor_pair, rnaseq or rnaseq_light) /!\ MANDATORY /!\."
  echo " -t <protocol>                    GenPipes protocol to be transfered (either: ensemble, sv or cancer). 'ensemble' and 'sv' are to be used with 'tumor_pair' <pipeline> and 'cancer' is to be used with 'rnaseq' <pipeline>"
  exit 1
  }

genpipes_tagging() {
  echo "Tagging GenPipes json"
  module load mugqic/python/3.11.1
  echo "$ABA_MOH/moh_automation/moh_automation_main/genpipes_deliverables_metrics.py -i $1 -o ${1/.json/_tagged.json}"
  module unload mugqic/python/3.11.1
}

genpipes_ingesting() {
  echo "Ingesting GenPipes json"
  source $ABA_MOH/project_tracking_cli/venv/bin/activate
  echo "pt-cli ingest genpipes --input-json $1"
  deactivate
}

genpipes_transfer() {
  echo "Transfering GenPipes"
  echo "$ABA_MOH/moh_automation/moh_automation_main/transfer_GenPipes.sh -r $1 -p $2 -t $3"
}

while getopts 'hj:r:p::t:' OPTION; do
  case "$OPTION" in
  # j)
  #   genpipes_json="$OPTARG"
  #   ;;
  r)
    readset_file="$OPTARG"
    # echo "readset_file: $readset_file"
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
  h)
    usage
    ;;
  ?)
    usage
    ;;
  esac
done

# mandatory arguments
if [ ! "$readset_file" ]; then
  echo -e "ERROR: Missing mandatory arguments -r.\n"
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

# Abacus MOH
ABA_MOH="/lb/project/mugqic/projects/MOH"
# ABA_MOH="/Users/pstretenowich/Mount_points/abacus/projects/MOH"
# Abacus main folder location
ABA_MAIN="$ABA_MOH/MAIN"

module load mugqic/genpipes
declare -A job_list_associative_array=()
{
  # Empty read to skip first line of readset file Cf. https://stackoverflow.com/questions/31911179/ignoring-first-line-column-header-while-reading-a-file-in-bash
  read -r
  while IFS=$'\t' read -r -a readset_fileArray; do
    sample=${readset_fileArray[0]}
    patient=$(echo "$sample" | sed -n 's/\(MoHQ-[CM|GC|HM|IQ|JG|MR|MU|XX]\{2\}-[0-9]\{1,\}-[a-zA-Z0-9]\{1,\}\).*/\1/p')
    if [ "$pipeline" = "tumor_pair" ]; then
      pipeline_job_list="TumorPair.${protocol}.job_list.*"
    elif [ "$pipeline" = "rnaseq" ]; then
      pipeline_job_list="RnaSeq.${protocol}.job_list.*"
    elif [ "$pipeline" = "rnaseq_light" ]; then
      pipeline_job_list="RnaSeqLight.job_list.*"
    fi
    # shellcheck disable=SC2086
    current_input_job_list=$(grep -m 1 -l "$patient" $ABA_MAIN/job_output/$pipeline_job_list | tr '\n' ' ')
    for file in $current_input_job_list; do
      if ! [[ -v job_list_associative_array["$file"] ]]; then
        job_list_associative_array["$file"]="$sample "
      else
        if ! [[ ${job_list_associative_array["$file"]} =~ (^| )$sample($| ) ]]; then
          job_list_associative_array["$file"]+="$sample "
        fi
      fi
    done
  done
} < "$readset_file"

input_jsons=""

for job_list in "${!job_list_associative_array[@]}"
do   
  # echo "key  : $job_list"
  # echo "value: ${job_list_associative_array[$job_list]}"
  failure=$(/Users/pstretenowich/Documents/local/apps/genpipes/utils/log_report.pl -nos $job_list | grep -v "^#")
  # echo "failure: $failure"
  if [[ -z $failure ]]; then
    # Let's tag GenPipes + Ingest GenPipes + Transfer it
    for sample in ${job_list_associative_array[$job_list]}; do
      # echo "sample: $sample"
      if [ "$pipeline" = "tumor_pair" ]; then
        pipeline_json="TumorPair.${protocol}_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].json"
      elif [ "$pipeline" = "rnaseq" ]; then
        pipeline_json="RnaSeq.${protocol}_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].json"
      elif [ "$pipeline" = "rnaseq_light" ]; then
        pipeline_json="RnaSeqLight_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].json"
      fi
      echo "pipeline_json: $pipeline_json"
      # shellcheck disable=SC2086
      input_jsons+="$(find $ABA_MAIN/json -type f -regex "$ABA_MAIN/json/$pipeline_json" -exec grep -m 1 -l "\"readset_name\": \"$sample" {} \; | tr '\n' ' ') "
      if [ -z "$input_jsons" ]; then
        echo "Nothing found in $ABA_MAIN/json/ for sample '$sample' and pipeline '$pipeline' and '$protocol'. Exiting."
        exit 0
      fi
      echo "input_jsons: $input_jsons"
    done
    for json in $input_jsons; do
      genpipes_tagging "$json"
      genpipes_ingesting "${json/.json/_tagged.json}"
      genpipes_transfer "$readset_file" "$pipeline" "$protocol"
    done
    echo "No failure found in $job_list"
  elif [[ $failure == *"FAILED"* ]]; then
    # Let's tag GenPipes + Ingest GenPipes
    for sample in ${job_list_associative_array[$job_list]}; do
      echo "sample: $sample"
      if [ "$pipeline" = "tumor_pair" ]; then
        pipeline_json="TumorPair.${protocol}_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].json"
      elif [ "$pipeline" = "rnaseq" ]; then
        pipeline_json="RnaSeq.${protocol}_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].json"
      elif [ "$pipeline" = "rnaseq_light" ]; then
        pipeline_json="RnaSeqLight_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].json"
      fi
      echo "pipeline_json: $pipeline_json"
      # shellcheck disable=SC2086
      input_jsons+="$(find $ABA_MAIN/json -type f -regex "$ABA_MAIN/json/$pipeline_json" -exec grep -m 1 -l "\"readset_name\": \"$sample" {} \; | tr '\n' ' ') "
      if [ -z "$input_jsons" ]; then
        echo "Nothing found in $ABA_MAIN/json/ for sample '$sample' and pipeline '$pipeline' and '$protocol'. Exiting."
        exit 0
      fi
      echo "input_jsons: $input_jsons"
    done
    for json in $input_jsons; do
      genpipes_tagging "$json"
      genpipes_ingesting "${json/.json/_tagged.json}"
    done
    echo "Failure found in $job_list"
  elif [[ $failure == *"ACTIVE"* ]]; then
    # Let's skip and wait
    echo "Job(s) still running for $job_list"
  else
    echo "Unknown status in $job_list"
  fi
done
