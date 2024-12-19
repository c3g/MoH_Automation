#!/usr/bin/env bash

THIS_SCRIPT=$(basename "$0")

usage() {
  echo "script usage: $THIS_SCRIPT -h [-l location] [-d destination] [-s sample]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -l <location>                    location of run to be transfered."
  echo " -d <destination>                 destination for the transfer (either Beluga or Cardinal or Abacus)."
  echo " -n <nucleic_acid_type>           nucleic_acid_type to be considered for the transfer (either DNA or RNA, Default: both)."
  echo " -s <sample>                      Path to file with Sample Name(s) (as they appear in the file <runid>-run.align_bwa_mem.csv) (by default it will consider ALL samples from the given run)."
  exit 1
  }

while getopts 'hl:d::n::s:' OPTION; do
  case "$OPTION" in
    l)
      location="$OPTARG"
      ;;
    d)
      destination="$OPTARG"
      ;;
    n)
      nucleic_acid_type="$OPTARG"
      ;;
    s)
      sample_file="$OPTARG"
      readarray -t sample < "$sample_file"
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
if [ ! "$location" ] || [ ! "$destination" ]; then
  echo -e "ERROR: Missing mandatory argument -l and/or -d.\n"
  usage
fi

if ! [[ $destination =~ Cardinal|Beluga|Abacus ]]; then
    echo -e "ERROR: Invalid destination: '$destination'. It has to be either Beluga, Cardinal or Abacus.\n"
    usage
fi

# location of processing data: Input to the script
if [[ $location == */ ]]
then
    echo "$location"
else
    location=$location"/"
fi

# # Getting client password to avoid getting timed out
# if [[ -n "${password}" ]]; then
#     echo -n Password: 
#     read -r -s password
#     echo
# fi

# The label is the run name based on the path given as argument
label=${location%?}
label=${label##*/}

# location on Beluga. CURRENTLY VERY IMPORTANT. DO NOT CHANGE OR IT WILL BREAK THE DATABASE
# SERIOUSLY DON'T CHANGE IT.
# PLEASE DONT.
if [[ $destination = Beluga ]]; then
    # Beluga main folder location
    DEST_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/raw_reads"
    # Beluga log file location
    DEST_LOG_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer"
    # Beluga run metrics location
    DEST_MET_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics"
    # Beluga Endpoint
    DEST_EP='278b9bfe-24da-11e9-9fa2-0a06afd4a22e'
elif [[ $destination = Cardinal ]]; then
    # Cardinal main folder location
    DEST_LOC="/project/60007/MOH/MAIN/raw_reads"
    # Cardinal log file location
    DEST_LOG_LOC="/project/60007/MOH/log_files/transfer"
    # Cardinal run metrics location
    DEST_MET_LOC="/project/60007/MOH/MAIN/metrics/run_metrics"
    # Cardinal Endpoint
    DEST_EP='26f926d9-6216-4e84-9037-a5c9567b5707'
elif [[ $destination = Abacus ]]; then
    # Abacus main folder location
    DEST_LOC="/lb/project/mugqic/projects/MOH/MAIN/raw_reads"
    # Abacus log file location
    DEST_LOG_LOC="/lb/project/mugqic/projects/MOH/log_files/transfer"
    # Abacus run metrics location
    DEST_MET_LOC="/lb/project/mugqic/projects/MOH/MAIN/metrics/run_metrics"
    # Abacus Endpoint
    DEST_EP=''
fi

# Temporary File location, you may want to change it to your scratch for easier clean up.
TEMP='/lb/project/mugqic/projects/MOH/TEMP'
TIMESTAMP=$(date +%FT%H.%M.%S)
LOGFILE="${label}_${TIMESTAMP}_transfer.log"
LISTFILE="${label}_${TIMESTAMP}_transfer.list"
touch "$TEMP/$LOGFILE"
touch "$TEMP/$LISTFILE"
echo "Log file of transfer from Abacus to $destination" > "$TEMP/$LOGFILE"

echo "-> Checking for files in Run $location"
echo "Transfered From $location" >> "$TEMP/$LOGFILE"

# abacus Endpoint
ABA_EP='26261fd6-0e6d-4252-a0ea-410b4b4f2eef'

# Loop over BAMS: DNA
if [ "$nucleic_acid_type" = "DNA" ] || [ -z "$nucleic_acid_type" ]; then
    if ls "$location"Aligned*/*/*/*/MoHQ*.bam 1> /dev/null 2>&1; then
        echo "Bams Transfered" >> "$TEMP/$LOGFILE"
        for i in "$location"Aligned*/*/*/*/MoHQ*.bam; do
            j=${i//.sorted.bam/.sorted.bai}
            sample_name=$(echo "$i" | cut -d'/' -f11)
            readset_name="$(echo "$i" | cut -d'/' -f11).$(echo "$i" | cut -d'/' -f12)"
            readset_name="${readset_name/run/}"
            file_name=$(echo "$i" | cut -d'/' -f13 | sed 's/.sorted.bam//g')
            LANE=$(echo "$i" | cut -d'/' -f9 | cut -d'.' -f2)
            DETAILS=$(awk -F "," -v LANE="$LANE" -v sample_name="$sample_name" '{ if ($3==LANE && $7==sample_name) {print $0}}' "$location"*run.csv)
            RUN_NAME=$(echo "$DETAILS" | awk -F "\"*,\"*" '{split($1, a, "_"); split(a[5], b, "-"); print b[1]}')
            RUNTYPE=$(echo "$DETAILS" | awk -F "\"*,\"*" '{print $4}')
            RUNID=$(echo "$DETAILS" | awk -F "\"*,\"*" '{print $2}')
            ADAP1=$(echo "$DETAILS" | awk -F "\"*,\"*" '{print $28}')
            ADAP2=$(echo "$DETAILS" | awk -F "\"*,\"*" '{print $29}')
            QUAL_OF=33
            BED=""
            BAM="raw_reads/$sample_name/${file_name}_${RUNID}_${LANE}.bam"
            # BAI="raw_reads/$sample_name/${file_name}_${RUNID}_${LANE}.bam.bai"
            FASTQ1=""
            FASTQ2=""
            if [[ ${#sample[@]} != 0 ]] && [[ ! " ${sample[*]} " =~ [[:space:]]${sample_name}[[:space:]] ]]; then
                :
            else
                # Make the oneliner readset file
                touch "$TEMP/${readset_name}_readset.tsv";
                echo -e 'Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM' > "$TEMP/${readset_name}_readset.tsv";
                echo -e "$sample_name\t${sample_name}.${RUNID}_${LANE}\t${RUNID}_${LANE}\t${RUNTYPE}\t${RUN_NAME}\t${LANE}\t${ADAP1}\t${ADAP2}\t${QUAL_OF}\t${BED}\t${FASTQ1}\t${FASTQ2}\t${BAM}" >> "$TEMP/${readset_name}_readset.tsv"
                {
                    # Adding readset to be transferred in a list file
                    echo "$TEMP/${readset_name}_readset.tsv $DEST_LOC/$sample_name/${readset_name}_readset.tsv"
                    # Adding bam and bai to be transferred in a list file
                    echo "$i $DEST_LOC/$sample_name/${file_name}_${RUNID}_${LANE}.bam"
                    echo "$j $DEST_LOC/$sample_name/${file_name}_${RUNID}_${LANE}.bam.bai"
                } >> "$TEMP/$LISTFILE"
                {
                    echo "$i,$sample_name/${file_name}_${RUNID}_${LANE}.bam"
                    echo "$j,$sample_name/${file_name}_${RUNID}_${LANE}.bam.bai"
                } >> "$TEMP/$LOGFILE"

            fi

        done;
    else
        echo "No BAM files found";
    fi;
fi;

# Loop over fastqs: DNA
if [ "$nucleic_acid_type" = "RNA" ] || [ -z "$nucleic_acid_type" ]; then
    if ls "$location"Unaligned*/*/*/MoHQ*_R1_001.fastq.gz 1> /dev/null 2>&1; then
        echo "fastq's Transfered" >> "$TEMP/$LOGFILE"
        for i in "$location"Unaligned*/*/*/MoHQ*_R1_001.fastq.gz; do
            j="${i/_R1_/_R2_}"
            sample_name=$(echo "$i" | cut -d'/' -f11 | cut -d'_' -f2)
            LANE=$(echo "$i" | cut -d'/' -f9 | cut -d'.' -f2)
            liba=$(echo "$i" | cut -d'/' -f8 | cut -d'_' -f2)
            libb=$(echo "$i" | cut -d'/' -f8 | cut -d'_' -f3)
            readset_name="${sample_name}.${liba}_${libb}_${LANE}"
            file_name1=$(echo "$i" | cut -d'/' -f12)
            file_name2="${file_name1/_R1_/_R2_}"
            DETAILS=$(awk -F "," -v LANE="$LANE" -v sample_name="$sample_name" '{ if ($3==LANE && $7==sample_name) {print $0}}' "$location"*run.csv)
            RUN_NAME=$(echo "$DETAILS" | awk -F "\"*,\"*" '{split($1, a, "_"); split(a[5], b, "-"); print b[1]}')
            RUNTYPE=$( echo "$DETAILS" | awk -F "\"*,\"*" '{print $4}' )
            RUNID=$( echo "$DETAILS" | awk -F "\"*,\"*" '{print $2}' )
            ADAP1=$( echo "$DETAILS" | awk -F "\"*,\"*" '{print $28}' )
            ADAP2=$( echo "$DETAILS" | awk -F "\"*,\"*" '{print $29}' )
            QUAL_OF=33
            BED=""
            BAM=""
            FASTQ1="raw_reads/$sample_name/$file_name1"
            FASTQ2="raw_reads/$sample_name/$file_name2"
            if [[ ${#sample[@]} != 0 ]] && [[ ! " ${sample[*]} " =~ [[:space:]]${sample_name}[[:space:]] ]]; then
                :
            else
                # Make the oneliner readset file
                touch "$TEMP/${readset_name}_readset.tsv"
                echo -e 'Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM' > "$TEMP/${readset_name}_readset.tsv"
                echo -e "$sample_name\t${sample_name}.${RUNID}_${LANE}\t${RUNID}_${LANE}\t${RUNTYPE}\t${RUN_NAME}\t${LANE}\t${ADAP1}\t${ADAP2}\t${QUAL_OF}\t${BED}\t${FASTQ1}\t${FASTQ2}\t${BAM}" >> "$TEMP/${readset_name}_readset.tsv"

                {
                    # Adding readset to be transferred in a list file
                    echo "$TEMP/${readset_name}_readset.tsv $DEST_LOC/$sample_name/${readset_name}_readset.tsv"
                    # Adding fastqs to be transferred in a list file
                    echo "$i $DEST_LOC/$sample_name/$file_name1"
                    echo "$j $DEST_LOC/$sample_name/$file_name2"
                } >> "$TEMP/$LISTFILE"
                {
                    echo "$i,$sample_name/$file_name1"
                    echo "$j,$sample_name/$file_name2"
                } >> "$TEMP/$LOGFILE"
            fi

        done;
    else
            echo "No fastq files found";
    fi;
fi;

echo "$TEMP/$LOGFILE $DEST_LOG_LOC/$LOGFILE" >> "$TEMP/$LISTFILE"

# Transfer over the metrics
MET_LOC=$(ls "$location"/*-novaseq*-run.align_bwa_mem.csv)
F_NAME=${MET_LOC##*/}

echo "$MET_LOC $DEST_MET_LOC/$F_NAME" >> "$TEMP/$LISTFILE"

if [[ $destination != Abacus ]]; then
    # Load globus module
    module load mugqic/globus-cli/3.24.0

    # Generate and store a UUID for the submission-id
    sub_id="$(globus task generate-submission-id)"

    # Start the batch transfer
    task_id="$(globus transfer --sync-level mtime --jmespath 'task_id' --format=UNIX --submission-id "$sub_id" --label "$label" --batch "$TEMP/$LISTFILE" $ABA_EP $DEST_EP)"

    echo "Waiting on 'globus transfer' task '$task_id'"
    globus task wait "$task_id" --polling-interval 60 -H
    # shellcheck disable=SC2181
    if [ $? -eq 0 ]; then
        module unload mugqic/globus-cli/3.24.0
        # shellcheck disable=SC1091
        source /lb/project/mugqic/projects/MOH/project_tracking_cli/venv/bin/activate
        # shellcheck disable=SC2086
        /lb/project/mugqic/projects/MOH/moh_automation/moh_automation_main/transfer2json.py --input $TEMP/$LISTFILE --source "abacus" --destination $destination --output /lb/project/mugqic/projects/MOH/Transfer_json/${LISTFILE/.list/.json} --operation_cmd_line "globus transfer --sync-level mtime --jmespath 'task_id' --format=UNIX --submission-id $sub_id --label $label --batch $TEMP/$LISTFILE $ABA_EP $DEST_EP"
        # if ! grep -q "password:" ~/.config/pt_cli/connect.yaml; then
        #     echo "  password: $password" >> ~/.config/pt_cli/connect.yaml
        # fi
        # shellcheck disable=SC2086
        pt-cli ingest transfer --input-json /lb/project/mugqic/projects/MOH/Transfer_json/${LISTFILE/.list/.json}
        # sed -i '/password: /d' ~/.config/pt_cli/connect.yaml
    else
        echo "$task_id failed!"
    fi
else
    while IFS= read -r line; do
        mkdir -p $(dirname $(echo $line | awk '{print $2}'))
        ln -sf $line
    done < "$TEMP/$LISTFILE"
fi
