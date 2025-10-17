#!/usr/bin/env bash

THIS_SCRIPT=$(basename "$0")

usage() {
  echo "script usage: $THIS_SCRIPT -h [-r readset_file] [-p pipeline] [-t protocol]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -r <readset_file>                Readset File used to submit GenPipes /!\ MANDATORY /!\."
  echo " -p <pipeline>                    GenPipes pipeline to be transfered (either: tumor_pair, rnaseq or rnaseq_light) /!\ MANDATORY /!\."
  echo " -t <protocol>                    GenPipes protocol to be transfered (either: ensemble, sv or cancer). 'ensemble' and 'sv' are to be used with 'tumor_pair' <pipeline> and 'cancer' is to be used with 'rnaseq' <pipeline>"
  exit 1
  }

while getopts 'hj:r:p::t:' OPTION; do
  case "$OPTION" in
  r)
    readset_file="$OPTARG"
    ;;
  p)
    pipeline="$OPTARG"
    case "$pipeline" in
      tumor_pair | rnaseq | rnaseq_light)
        ;;
      *)
        echo -e "ERROR: Invalid pipeline: '$pipeline'.\n"
        usage
        ;;
    esac
    ;;
  t)
    protocol="$OPTARG"
    if [ -n "$protocol" ]; then
        case "$protocol" in
          ensemble | sv | cancer)
            ;;
          *)
            echo -e "ERROR: Invalid protocol: '$protocol'.\n"
            usage
            ;;
        esac
      fi
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

# The label is the readset file name
label=${readset_file%.*}
label=${label##*/}
label="${label}.${pipeline}.${protocol}"

# Destination cluster
## Beluga main folder location
# DEST_MAIN="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN"
## Beluga log file location
# DEST_LOG_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer"
## Beluga Endpoint
# DEST_EP='278b9bfe-24da-11e9-9fa2-0a06afd4a22e'
## Rorqual main folder location
DEST_MAIN="/project/6007512/C3G/projects/MOH_PROCESSING/MAIN"
## Rorqual Endpoint
DEST_EP='f19f13f5-5553-40e3-ba30-6c151b9d35d4'

# Source cluster
if [[ $HOSTNAME == "abacus"* ]]; then
  # Abacus MOH
  SRC_MOH="/lb/project/mugqic/projects/MOH"
  # Abacus main folder location
  SRC_MAIN="$SRC_MOH/MAIN"
  # Abacus log file location
  SRC_LOG_LOC="$SRC_MOH/TEMP"
  # Abacus Endpoint
  SRC_EP='26261fd6-0e6d-4252-a0ea-410b4b4f2eef'
  # Abacus name for json creation
  SRC='abacus'
elif [[ $HOSTNAME == "cardinal"* ]]; then
  # Cardinal MOH
  SRC_MOH="/project/def-c3g/MOH"
  # Cardinal main folder location
  SRC_MAIN="$SRC_MOH/MAIN"
  # Cardinal log file location
  SRC_LOG_LOC="$SRC_MOH/log_files/transfer"
  # Cardinal Endpoint
  SRC_EP='a6df16bc-4e7f-4784-9afa-8ceb7b20b7c0'
  # Cardinal name for json creation
  SRC='cardinal'
else
  echo "ERROR: Unknown cluster. Exiting."
  exit 1
fi

TIMESTAMP=$(date +%FT%H.%M.%S)
# LOGFILE="$TIMESTAMP_transfer_GenPipes.log"
LISTFILE="${label}_${TIMESTAMP}_transfer_GenPipes.list"
merged_json="${label}_${TIMESTAMP}_transfer_GenPipes_files.json"
timestamp_start=$(date "+%Y-%m-%dT%H.%M.%S")

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
    patient=$(echo "$sample" | sed -n 's/\(MoHQ-[CM|GC|HM|IQ|JG|MR|MU|XX]\{2\}-[0-9]\{1,\}-[a-zA-Z0-9]\{1,\}\).*/\1/p')
    patients_associative_array["$patient"]+="$sample "
    samples_associative_array["$sample"]+="$readset "
    if [ "$pipeline" = "tumor_pair" ]; then
      pipeline_json="TumorPair.${protocol}_*_tagged.json"
    elif [ "$pipeline" = "rnaseq" ]; then
      pipeline_json="RnaSeq.${protocol}_*_tagged.json"
    elif [ "$pipeline" = "rnaseq_light" ]; then
      pipeline_json="RnaSeqLight_*_tagged.json"
    fi
    # shellcheck disable=SC2086
    input_jsons=$(grep -m 1 -l "\"readset_name\": \"$sample" $SRC_MAIN/json/$pipeline_json | tr '\n' ' ')
    if [ -z "$input_jsons" ]; then
      echo "Nothing found in $SRC_MAIN/json/ for sample '$sample' and pipeline '$pipeline' and protocol '$protocol'. Exiting."
      exit 0
    fi
    # shellcheck disable=SC2086
    $SRC_MOH/moh_automation/moh_automation_main/merge_GenPipes_jsons.py -i $input_jsons -o $SRC_MOH/Transfer_json/$merged_json
    echo "Merged GenPipes JSONs written to $SRC_MOH/Transfer_json/$merged_json"
  done
} < "$readset_file"
module unload mugqic/python/3.10.4

for patient in "${!patients_associative_array[@]}"; do
  if [ "$pipeline" = "tumor_pair" ] && [ "$protocol" = "ensemble" ]; then
    {
      echo "--recursive $SRC_MAIN/pairedVariants/${patient} $DEST_MAIN/pairedVariants/${patient}"
      echo "--recursive $SRC_MAIN/pairedVariants/ensemble/${patient} $DEST_MAIN/pairedVariants/ensemble/${patient}"
      echo "--recursive $SRC_MAIN/SVariants/${patient} $DEST_MAIN/SVariants/${patient}"
      echo "$SRC_MAIN/metrics/dna/${patient}.multiqc.html $DEST_MAIN/metrics/dna/${patient}.multiqc.html"
      echo "--recursive $SRC_MAIN/metrics/dna/${patient}.multiqc_data $DEST_MAIN/metrics/dna/${patient}.multiqc_data"

    } >> "$SRC_LOG_LOC/$LISTFILE"
  elif [ "$pipeline" = "tumor_pair" ] && [ "$protocol" = "sv" ]; then
    {
      echo "--recursive $SRC_MAIN/SVariants/${patient} $DEST_MAIN/SVariants/${patient}"

    } >> "$SRC_LOG_LOC/$LISTFILE"
  fi


  for sample in ${patients_associative_array[$patient]} ; do
    if [ "$pipeline" = "rnaseq_light" ]; then
      {
        echo "--recursive $SRC_MAIN/kallisto/${sample} $DEST_MAIN/kallisto/${sample}"
        echo "--recursive $SRC_MAIN/trim/${sample} $DEST_MAIN/trim/${sample}"

      } >> "$SRC_LOG_LOC/$LISTFILE"
    fi

    if [ "$pipeline" = "rnaseq" ] && [ "$protocol" = "cancer" ]; then
      {
        echo "--recursive $SRC_MAIN/alignment/${sample} $DEST_MAIN/alignment/${sample}"
        echo "--recursive $SRC_MAIN/fusion/${sample} $DEST_MAIN/fusion/${sample}"
        echo "--recursive $SRC_MAIN/metrics/${sample} $DEST_MAIN/metrics/${sample}"
        echo "--recursive $SRC_MAIN/metrics/multiqc_by_sample/${sample} $DEST_MAIN/metrics/multiqc_by_sample/${sample}"
        echo "--recursive $SRC_MAIN/metrics/sortmerna/${sample} $DEST_MAIN/metrics/sortmerna/${sample}"
        echo "$SRC_MAIN/tracks/bigWig/${sample}.bw $DEST_MAIN/tracks/bigWig/${sample}.bw"
        echo "--recursive $SRC_MAIN/trim/${sample} $DEST_MAIN/trim/${sample}"

      } >> "$SRC_LOG_LOC/$LISTFILE"
    fi

    if [ "$pipeline" = "tumor_pair" ] && [ "$protocol" = "ensemble" ]; then
      {
        echo "--recursive $SRC_MAIN/alignment/${sample} $DEST_MAIN/alignment/${sample}"
        echo "--recursive $SRC_MAIN/metrics/dna/${sample} $DEST_MAIN/metrics/dna/${sample}"

      } >> "$SRC_LOG_LOC/$LISTFILE"
      if [[ $sample == *DT ]]; then
        {
          echo "$SRC_MAIN/metrics/${sample}.concordance.tsv $DEST_MAIN/metrics/${sample}.concordance.tsv"
          echo "$SRC_MAIN/metrics/${sample}.contamination.tsv $DEST_MAIN/metrics/${sample}.contamination.tsv"

        } >> "$SRC_LOG_LOC/$LISTFILE"
      fi
    fi
  done
done

echo "List of files to transfer written to $SRC_LOG_LOC/$LISTFILE"

# Load globus module
module load mugqic/globus-cli/3.24.0

# Generate and store a UUID for the submission-id
sub_id="$(globus task generate-submission-id)"

# Start the batch transfer
task_id="$(globus transfer --sync-level mtime --jmespath 'task_id' --format=UNIX --submission-id "$sub_id" --label "$label" --batch "$SRC_LOG_LOC/$LISTFILE" $SRC_EP $DEST_EP)"

echo "Waiting on 'globus transfer' task '$task_id'"
globus task wait "$task_id" --polling-interval 60 -H
# shellcheck disable=SC2181
if [ $? -eq 0 ]; then
  transfer_json=$SRC_MOH/Transfer_json/${LISTFILE/.list/.json}
  module unload mugqic/globus-cli/3.24.0
  timestamp_end=$(date "+%Y-%m-%dT%H.%M.%S")
  # shellcheck disable=SC1091
  source $SRC_MOH/project_tracking_cli/venv/bin/activate
  # shellcheck disable=SC2086
  $SRC_MOH/moh_automation/moh_automation_main/transfer2json.py --input $SRC_LOG_LOC/$LISTFILE --output $transfer_json --source $SRC --destination rorqual --operation_cmd_line "globus transfer --sync-level mtime --jmespath 'task_id' --format=UNIX --submission-id $sub_id --label $label --batch $SRC_LOG_LOC/$LISTFILE $SRC_EP $DEST_EP" --genpipes $SRC_MOH/Transfer_json/$merged_json --start $timestamp_start --stop $timestamp_end
  # shellcheck disable=SC2086
  pt-cli ingest transfer --input-json $transfer_json
else
    echo "$task_id failed!"
fi
