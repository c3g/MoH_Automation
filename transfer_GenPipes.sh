#!/usr/bin/env bash

usage() {
  echo "script usage: transfer_GenPipes.sh -h [-j genpipes_json] [-r readset_file]"
  echo "Usage:"
  echo " -h                               Display this help message."
  # echo " -j <genpipes_json>               json file generated by GenPipes /!\ MANDATORY /!\."
  echo " -r <readset_file>                Readset File used to submit GenPipes /!\ MANDATORY /!\."
  echo " -p <pipeline>                    GenPipes pipeline to be transfered (either: tumor_pair, rnaseq or rnaseq_light) /!\ MANDATORY /!\."
  echo " -t <protocol>                    GenPipes protocol to be transfered (either: ensemble, sv or cancer). 'ensemble' and 'sv' are to be used with 'tumor_pair' <pipeline> and 'cancer' is to be used with 'rnaseq' <pipeline>"
  exit 1
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

# operation_cmd_line=$(jq '.operation_cmd_line' "$genpipes_json")
# pipeline=$(echo "$operation_cmd_line" | cut -d' ' -f1 | rev | cut -d'/' -f2 | rev)
# protocol=$(echo "$genpipes_json" | cut -d'.' -f2 |  cut -d'_' -f1)

# mandatory arguments
if [ ! "$readset_file" ]; then
  echo -e "ERROR: Missing mandatory arguments -r.\n"
  usage
fi
# if [ ! "$readset_file" ] || [ ! "$genpipes_json" ]; then
#   echo -e "ERROR: Missing mandatory arguments -r and -j.\n"
#   usage
# fi

# wrong protocol with wrong pipeline
if [ "$pipeline" = "tumor_pair" ] && ! [[ "$protocol" == "ensemble" || "$protocol" == "sv" ]]; then
  echo -e "ERROR: pipeline: '$pipeline' only accepts protocol: 'ensemble' or 'sv', '$protocol' provided.\n"
  usage
elif [ "$pipeline" = "rnaseq" ] && [ "$protocol" != "cancer" ]; then
  echo -e "ERROR: pipeline: '$pipeline' only accepts protocol: 'cancer', '$protocol' provided.\n"
  usage
fi

# The label is the readset file name
label=${readset_file%.*}
label=${label##*/}

# Beluga main folder location
BEL_MAIN="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN"
# Beluga log file location
# BEL_LOG_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer"
# Beluga Endpoint
BEL_EP='278b9bfe-24da-11e9-9fa2-0a06afd4a22e'
# Abacus MOH
ABA_MOH="/lb/project/mugqic/projects/MOH"
# Abacus main folder location
ABA_MAIN="$ABA_MOH/MAIN"
# Abacus log file location
ABA_LOG_LOC="$ABA_MOH/TEMP"
# Abacus Endpoint
ABA_EP='26261fd6-0e6d-4252-a0ea-410b4b4f2eef'

TIMESTAMP=$(date +%FT%H.%M.%S)
# LOGFILE="$TIMESTAMP_transfer_GenPipes.log"
LISTFILE="${label}_${TIMESTAMP}_transfer_GenPipes.list"
merged_json="${label}_${TIMESTAMP}_transfer_GenPipes_files.json"

# samples files, case,normal,tumor
# Activate beluga endpoint, before

declare -A patients_associative_array=()
declare -A samples_associative_array=()

module load mugqic/python/3.10.4
{
  # Empty read to skip first line of readset file Cf. https://stackoverflow.com/questions/31911179/ignoring-first-line-column-header-while-reading-a-file-in-bash
  read -r
  while IFS=$'\t' read -r -a readset_fileArray; do
    sample=${readset_fileArray[0]}
    readset=${readset_fileArray[1]}
    patient=$(echo "$sample" | sed -n 's/\(MoHQ-[JG|HM|CM|GC|MU|MR|XX]\{2\}-[0-9]\{1,\}-[0-9]\{1,\}\).*/\1/p')
    patients_associative_array["$patient"]+="$sample "
    samples_associative_array["$sample"]+="$readset "
    if [ "$pipeline" = "tumor_pair" ]; then
      pipeline_json="TumorPair.${protocol}_*.json"
    elif [ "$pipeline" = "rnaseq" ]; then
      pipeline_json="RnaSeq.${protocol}_*.json"
    elif [ "$pipeline" = "rnaseq_light" ]; then
      pipeline_json="RnaSeqLight_*.json"
    fi
    # shellcheck disable=SC2086
    input_jsons=$(grep -l "\"readset_name\": \"$sample" $ABA_MAIN/json/$pipeline_json | tr '\n' ' ')
    if [ -n "$input_jsons" ]; then
      echo "Nothing found in $ABA_MAIN/json/ for sample '$sample' and pipeline '$pipeline' and '$protocol'. Exiting."
      exit 0
    fi
    # shellcheck disable=SC2086
    $ABA_MOH/moh_automation/moh_automation_main/merge_GenPipes_jsons.py -i $input_jsons -o $ABA_MOH/Transfer_json/$merged_json
  done
} < "$readset_file"
module unload mugqic/python/3.10.4

for patient in "${!patients_associative_array[@]}"; do
  if [ "$pipeline" = "tumor_pair" ] && [ "$protocol" = "ensemble" ]; then
    {
      echo "--recursive $ABA_MAIN/alignment/realign/${patient} $BEL_MAIN/alignment/realign/${patient}"
      echo "--recursive $ABA_MAIN/pairedVariants/${patient} $BEL_MAIN/pairedVariants/${patient}"
      echo "--recursive $ABA_MAIN/SVariants/${patient} $BEL_MAIN/SVariants/${patient}"
      echo "$ABA_MAIN/metrics/dna/${patient}.multiqc.html $BEL_MAIN/metrics/dna/${patient}.multiqc.html"
      echo "--recursive $ABA_MAIN/metrics/dna/${patient}.multiqc_data $BEL_MAIN/metrics/dna/${patient}.multiqc_data"

    } >> "$ABA_LOG_LOC/$LISTFILE"
  elif [ "$pipeline" = "tumor_pair" ] && [ "$protocol" = "sv" ]; then
    {
      echo "--recursive $ABA_MAIN/SVariants/${patient} $BEL_MAIN/SVariants/${patient}"

    } >> "$ABA_LOG_LOC/$LISTFILE"
  fi


  for sample in ${patients_associative_array[$patient]} ; do
    if [ "$pipeline" = "rnaseq_light" ]; then
      {
        echo "--recursive $ABA_MAIN/kallisto/${sample} $BEL_MAIN/kallisto/${sample}"
        echo "--recursive $ABA_MAIN/trim/${sample} $BEL_MAIN/trim/${sample}"

      } >> "$ABA_LOG_LOC/$LISTFILE"
    fi

    if [ "$pipeline" = "rnaseq" ] && [ "$protocol" = "cancer" ]; then
      {
        echo "--recursive $ABA_MAIN/alignment/${sample} $BEL_MAIN/alignment/${sample}"
        echo "--recursive $ABA_MAIN/fusion/${sample} $BEL_MAIN/fusion/${sample}"
        echo "--recursive $ABA_MAIN/metrics/${sample} $BEL_MAIN/metrics/${sample}"
        echo "--recursive $ABA_MAIN/metrics/multiqc_by_sample/${sample} $BEL_MAIN/metrics/multiqc_by_sample/${sample}"
        echo "--recursive $ABA_MAIN/metrics/sortmerna/${sample} $BEL_MAIN/metrics/sortmerna/${sample}"
        echo "$ABA_MAIN/tracks/bigWig/${sample}.bw $BEL_MAIN/tracks/bigWig/${sample}.bw"
        echo "--recursive $ABA_MAIN/trim/${sample} $BEL_MAIN/trim/${sample}"

      } >> "$ABA_LOG_LOC/$LISTFILE"
    fi

    if [ "$pipeline" = "tumor_pair" ] && [ "$protocol" = "ensemble" ]; then
      {
        echo "--recursive $ABA_MAIN/alignment/${sample} $BEL_MAIN/alignment/${sample}"
        echo "--recursive $ABA_MAIN/metrics/dna/${sample} $BEL_MAIN/metrics/dna/${sample}"

      } >> "$ABA_LOG_LOC/$LISTFILE"
      if [[ $sample == *DT ]]; then
        {
          echo "$ABA_MAIN/metrics/${sample}.concordance.tsv $BEL_MAIN/metrics/${sample}.concordance.tsv"
          echo "$ABA_MAIN/metrics/${sample}.contamination.tsv $BEL_MAIN/metrics/${sample}.contamination.tsv"

        } >> "$ABA_LOG_LOC/$LISTFILE"
      fi
    fi
  done
done

# Load globus module
module load mugqic/globus-cli/3.24.0

# Generate and store a UUID for the submission-id
sub_id="$(globus task generate-submission-id)"

# Start the batch transfer
task_id="$(globus transfer --jmespath 'task_id' --format=UNIX --submission-id "$sub_id" --label "$label" --batch "$ABA_LOG_LOC/$LISTFILE" $ABA_EP $BEL_EP)"

echo "Waiting on 'globus transfer' task '$task_id'"
globus task wait "$task_id" --polling-interval 60 -H
# shellcheck disable=SC2181
if [ $? -eq 0 ]; then
  module unload mugqic/globus-cli/3.24.0
  # shellcheck disable=SC1091
  source $ABA_MOH/project_tracking_cli/venv/bin/activate
  # shellcheck disable=SC2086
  $ABA_MOH/moh_automation/moh_automation_main/transfer2json.py --input $ABA_LOG_LOC/$LISTFILE --output $ABA_MOH/Transfer_json/${LISTFILE/.txt/.json} --operation_cmd_line "globus transfer --submission-id $sub_id --label $label --batch $ABA_LOG_LOC/$LISTFILE $ABA_EP $BEL_EP" --genpipes $ABA_MOH/Transfer_json/$merged_json
  # shellcheck disable=SC2086
  pt-cli ingest transfer --input-json $ABA_MOH/Transfer_json/${LISTFILE/.txt/.json}
else
    echo "$task_id failed!"
fi
