#!/usr/bin/env bash

usage() {
  echo "script usage: transfer_GenPipes.sh -h [-r readset_file] [-p pipeline] [-t protocol]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -r <readset_file>                Readset File used to submit GenPipes /!\ MANDATORY /!\."
  echo " -p <pipeline>                    GenPipes pipeline to be transfered (either: tumor_pair, rnaseq or rnaseq_light) /!\ MANDATORY /!\."
  echo " -t <protocol>                    GenPipes protocol to be transfered (either: ensemble, sv or cancer). 'ensemble' and 'sv' are to be used with 'tumor_pair' <pipeline> and 'cancer' is to be used with 'rnaseq' <pipeline>"
  exit 1
  }

while getopts 'hr:p::t:' OPTION; do
  case "$OPTION" in
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
if [ ! "$readset_file" ] || [ ! "$pipeline" ]; then
  echo -e "ERROR: Missing mandatory arguments -r and/or -p.\n"
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

# Beluga main folder location
BEL_MAIN="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN"
# Beluga log file location
# BEL_LOG_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer"
# Beluga Endpoint
BEL_EP='278b9bfe-24da-11e9-9fa2-0a06afd4a22e'
# Abacus main folder location
ABA_MAIN="/lb/project/mugqic/projects/MOH/MAIN"
# Abacus log file location
ABA_LOG_LOC="/lb/project/mugqic/projects/MOH/TEMP"
# Abacus Endpoint
ABA_EP='26261fd6-0e6d-4252-a0ea-410b4b4f2eef'

TIMESTAMP=$(date +%FT%H.%M.%S)
# LOGFILE="$TIMESTAMP_transfer_GenPipes.log"
LISTFILE="${TIMESTAMP}_transfer_GenPipes.list"

# samples files, case,normal,tumor
# Activate beluga endpoint, before
readset_file=$1
pipeline=$1

declare -A patients_associative_array=()
declare -A samples_associative_array=()

{
  while IFS=$'\t' read -r -a readset_fileArray; do
    sample=${readset_fileArray[0]}
    readset=${readset_fileArray[1]}
    patient=$(echo "$sample" | sed -n 's/\(MoHQ-[JG|HM|CM|GC|MU|MR|XX]\{2\}-[0-9]\{1,\}-[0-9]\{1,\}\).*/\1/p')
    patients_associative_array["$patient"]+="$sample "
    samples_associative_array["$sample"]+="$readset "
  done
} < "$readset_file"

for patient in "${!patients_associative_array[@]}"; do
  echo "$patient"
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
    echo "$sample"
    if [ "$pipeline" = "rnaseq_light" ]; then
      {
        echo "--recursive $ABA_MAIN/kallisto/${sample} $BEL_MAIN/kallisto/${sample}"
        echo "--recursive $ABA_MAIN/trim/${sample} $BEL_MAIN/trim/${sample}"

      } >> "$ABA_LOG_LOC/$LISTFILE"
    fi

    if [ "$pipeline" = "rnaseq" ] && [ "$protocol" = "cancer" ]; then
      {
        echo "--recursive $ABA_MAIN/alignment/${sample} $BEL_MAIN/alignment/${sample}"
        echo "--recursive $ABA_MAIN/alignment_1stPass/${sample} $BEL_MAIN/alignment_1stPass/${sample}"
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
label=${readset_file%.*}
label=${label##*/}

# Start the batch transfer
task_id="$(globus transfer --jmespath 'task_id' --format=UNIX --submission-id "$sub_id" --label "$label" --batch "$TEMP/$LISTFILE" $ABA_EP $BEL_EP)"

echo "Waiting on 'globus transfer' task '$task_id'"
globus task wait "$task_id" --polling-interval 60 -H
# shellcheck disable=SC2181
if [ $? -eq 0 ]; then
  module unload mugqic/globus-cli/3.24.0
  # shellcheck disable=SC1091
  source /lb/project/mugqic/projects/MOH/project_tracking_cli/venv/bin/activate
  # shellcheck disable=SC2086
  /lb/project/mugqic/projects/MOH/moh_automation/moh_automation_main/transfer2json.py --input $TEMP/$LISTFILE --output /lb/project/mugqic/projects/MOH/Transfer_json/${LISTFILE/.txt/.json} --operation_cmd_line "globus transfer --submission-id $sub_id --label $label --batch $TEMP/$LISTFILE $ABA_EP $BEL_EP"
  # shellcheck disable=SC2086
  pt-cli ingest transfer --input-json /lb/project/mugqic/projects/MOH/Transfer_json/${LISTFILE/.txt/.json}
else
    echo "$task_id failed!"
fi
