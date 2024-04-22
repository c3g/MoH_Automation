#!/usr/bin/env bash

THIS_SCRIPT=$(basename "$0")

usage() {
  echo "script usage: $THIS_SCRIPT -h [-l location] [-s sample]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -l <location>                    location of run to be transfered."
  echo " -s <sample>                      Path to file with Sample Name(s) (as they appear in the file <runid>-run.align_bwa_mem.csv) (by default it will consider ALL samples from the given run)."
  exit 1
  }

while getopts 'hl::s:' OPTION; do
  case "$OPTION" in
    l)
      location="$OPTARG"
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
if [ ! "$location" ]; then
  echo -e "ERROR: Missing mandatory arguments -l.\n"
  usage
fi

# location of processing data: Input to the script
if [[ $location == */ ]]
then
    echo "$location"
else
    location=$location"/"
fi

# Getting client password to avoid getting timed out
echo -n Password: 
read -r -s password

# The label is the run name based on the path given as argument
label=${location%?}
label=${label##*/}

# Temporary File location, you may want to change it to your scratch for easier clean up.
TEMP='/lb/project/mugqic/projects/MOH/TEMP'
TIMESTAMP=$(date +%FT%H.%M.%S)
LOGFILE="${label}_${TIMESTAMP}_transfer.log"
LISTFILE="${label}_${TIMESTAMP}_transfer.list"
touch "$TEMP/$LOGFILE"
touch "$TEMP/$LISTFILE"
echo "Log file of transfer from Abacus to Beluga" > "$TEMP/$LOGFILE"

echo "-> Checking for files in Run $location"
echo "Transfered From $location" >> "$TEMP/$LOGFILE"

# location on Beluga. CURRENTLY VERY IMPORTANT. DO NOT CHANGE OR IT WILL BREAK THE DATABASE
# SERIOUSLY DON'T CHANGE IT.
# PLEASE DONT.
BEL_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/raw_reads/"
# Beluga log file location
BEL_LOG_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/"
# Beluga Endpoint
BEL_EP='278b9bfe-24da-11e9-9fa2-0a06afd4a22e'
# Narval Endpoint #not used
# NAR_EP='a1713da6-098f-40e6-b3aa-034efe8b6e5b'
# abacus Endpoint
ABA_EP='26261fd6-0e6d-4252-a0ea-410b4b4f2eef'

# Loop over BAMS
# This is specific for MoHQ named Bams!
if ls "$location"Aligned*/*/*/*/MoHQ*.bam 1> /dev/null 2>&1; then
    echo "Bams Transfered" >> "$TEMP/$LOGFILE"
    for i in "$location"Aligned*/*/*/*/MoHQ*.bam; do
        j=${i//.sorted.bam/.sorted.bai}
        # j=$(echo "$i" | sed 's/.sorted.bam/.sorted.bai/g')
        sample_name=$(echo "$i" | cut -d'/' -f11)
        readset_name="$(echo "$i" | cut -d'/' -f11).$(echo "$i" | cut -d'/' -f12)"
        readset_name="${readset_name/run/}"
        file_name=$(echo "$i" | cut -d'/' -f13 | sed 's/.sorted.bam//g')
        DETAILS=$(cat "$location"*run.csv | grep "$sample_name")
        RUN_NAME=$(echo "$DETAILS" | awk -F "\"*,\"*" '{split($1, a, "_"); split(a[5], b, "-"); print b[1]}')
        RUNTYPE=$(echo "$DETAILS" | awk -F "\"*,\"*" '{print $4}')
        RUNID=$(echo "$DETAILS" | awk -F "\"*,\"*" '{print $2}')
        LANE=$(echo "$DETAILS" | awk -F "\"*,\"*" '{print $3}')
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
                echo "$TEMP/${readset_name}_readset.tsv $BEL_LOC$sample_name/${readset_name}_readset.tsv"
                # Adding bam and bai to be transferred in a list file
                echo "$i $BEL_LOC$sample_name/${file_name}_${RUNID}_${LANE}.bam"
                echo "$j $BEL_LOC$sample_name/${file_name}_${RUNID}_${LANE}.bam.bai"
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

# RNA now or possibly both
if ls "$location"Unaligned*/*/*/MoHQ*_R1_001.fastq.gz 1> /dev/null 2>&1; then
    echo "fastq's Transfered" >> "$TEMP/$LOGFILE"
    for i in "$location"Unaligned*/*/*/MoHQ*_R1_001.fastq.gz; do
        j="${i/_R1_/_R2_}"
        sample_name=$(echo "$i" | cut -d'/' -f11 | cut -d'_' -f2)
        lane=$(echo "$i" | cut -d'/' -f9 | cut -d'.' -f2)
        liba=$(echo "$i" | cut -d'/' -f8 | cut -d'_' -f2)
        libb=$(echo "$i" | cut -d'/' -f8 | cut -d'_' -f3)
        readset_name="${sample_name}.${liba}_${libb}_${lane}"
        file_name1=$(echo "$i" | cut -d'/' -f12)
        file_name2="${file_name1/_R1_/_R2_}"
        DETAILS=$(cat "$location"*run.csv | grep "$sample_name")
        RUN_NAME=$(echo "$DETAILS" | awk -F "\"*,\"*" '{split($1, a, "_"); split(a[5], b, "-"); print b[1]}')
        RUNTYPE=$( echo "$DETAILS" | awk -F "\"*,\"*" '{print $4}' )
        RUNID=$( echo "$DETAILS" | awk -F "\"*,\"*" '{print $2}' )
        LANE=$( echo "$DETAILS" | awk -F "\"*,\"*" '{print $3}' )
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
                echo "$TEMP/${readset_name}_readset.tsv $BEL_LOC$sample_name/${readset_name}_readset.tsv"
                # Adding fastqs to be transferred in a list file
                echo "$i $BEL_LOC$sample_name/$file_name1"
                echo "$j $BEL_LOC$sample_name/$file_name2"
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

echo "$TEMP/$LOGFILE $BEL_LOG_LOC$LOGFILE" >> "$TEMP/$LISTFILE"

# Transfer over the metrics
MET_LOC=$(ls "$location"/*-novaseq-run.align_bwa_mem.csv)
F_NAME=${MET_LOC##*/}

echo "$MET_LOC /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/$F_NAME" >> "$TEMP/$LISTFILE"

# Load globus module
module load mugqic/globus-cli/3.24.0

# Generate and store a UUID for the submission-id
sub_id="$(globus task generate-submission-id)"

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
    echo "  password: $password" >> ~/.config/pt_cli/connect.yaml
    # shellcheck disable=SC2086
    pt-cli ingest transfer --input-json /lb/project/mugqic/projects/MOH/Transfer_json/${LISTFILE/.txt/.json}
    sed -i '/password: /d' ~/.config/pt_cli/connect.yaml
else
    echo "$task_id failed!"
fi