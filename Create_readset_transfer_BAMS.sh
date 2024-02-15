#!/usr/bin/env bash

# Temporary File Location, you may want to change it to your scratch for easier clean up.
TEMP='/lb/project/mugqic/projects/MOH/TEMP'
TIMESTAMP=$(date +%FT%H.%M.%S)
LOGFILE=$TIMESTAMP"_transfer_log.txt"
LISTFILE=$TIMESTAMP"_transfer_list.txt"
touch "$TEMP/$LOGFILE"
touch "$TEMP/$LISTFILE"
echo "Log file of transfer from abacus to Beluga" > "$TEMP/$LOGFILE"

# Location of processing data
Location="$1"
if [[ $Location == */ ]]
then
    echo "$Location"
else
    Location=$Location"/"
fi
echo "$Location"
echo "Transfered From $Location" >> "$TEMP/$LOGFILE"

# Location on Beluga. CURRENTLY VERY IMPORTANT. DO NOT CHANGE OR IT WILL BREAK THE DATABASE
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
if ls "$Location"Aligned*/*/*/*/MoHQ*.bam 1> /dev/null 2>&1; then
    echo "Bams Transfered" >> "$TEMP/$LOGFILE"
    for i in "$Location"Aligned*/*/*/*/MoHQ*.bam; do
        j=${i//.sorted.bam/.sorted.bai}
        # j=$(echo "$i" | sed 's/.sorted.bam/.sorted.bai/g')
    	sample_name=$(echo "$i" | cut -d'/' -f11)
        readset_name="$(echo "$i" | cut -d'/' -f11)_$(echo "$i" | cut -d'/' -f12)"
        file_name=$(echo "$i" | cut -d'/' -f13 | sed 's/.sorted.bam//g')
        # Make the oneliner readset file
        touch "$TEMP/${readset_name}_readset.tsv";
        echo -e 'Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM' > "$TEMP/${readset_name}_readset.tsv";
        DETAILS=$(cat "$Location"*run.csv | grep "$sample_name")
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

    done;
else
    echo "No BAM files found";
fi;

# RNA now or possibly both
if ls "$Location"Unaligned*/*/*/MoHQ*_R1_001.fastq.gz 1> /dev/null 2>&1; then
    echo "fastq's Transfered" >> "$TEMP/$LOGFILE"
    for i in "$Location"Unaligned*/*/*/MoHQ*_R1_001.fastq.gz; do
        j="${i/_R1_/_R2_}"
        sample_name=$(echo "$i" | cut -d'/' -f11)
        sample_name="${sample_name%_*}"
        sample_name="${sample_name/Sample_/}"
        readset_name=$(echo "$i" | cut -d'/' -f11)
        readset_name="${readset_name/Sample_/}"
        file_name1=$(echo "$i" | cut -d'/' -f12)
        file_name2="${file_name1/_R1_/_R2_}"
        # Make the oneliner readset file
        touch "$TEMP/${readset_name}_readset.tsv"
        echo -e 'Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM' > "$TEMP/${readset_name}_readset.tsv"
        DETAILS=$(cat "$Location"*run.csv | grep "$sample_name")
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

    done;
else
        echo "No fastq files found";
fi;

echo "$TEMP/$LOGFILE $BEL_LOG_LOC$LOGFILE" >> "$TEMP/$LISTFILE"

# Transfer over the metrics
MET_LOC=$(ls "$Location"/*-novaseq-run.align_bwa_mem.csv)
F_NAME=${MET_LOC##*/}

echo "$MET_LOC /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/$F_NAME" >> "$TEMP/$LISTFILE"

# Load globus module
module load mugqic/globus-cli/3.7.0

# Generate and store a UUID for the submission-id
sub_id="$(globus task generate-submission-id)"
label=${Location%?}
label=${label##*/}

# Start the batch transfer
task_id="$(globus transfer --submission-id "$sub_id" --label "$label" --batch "$TEMP/$LISTFILE" $ABA_EP $BEL_EP)"

echo "Waiting on 'globus transfer' task '$task_id'"
globus task wait "$task_id" --timeout 30
if [ $? -eq 0 ]; then
    echo "source /lb/project/mugqic/projects/MOH/project_tracking_cli/venv/bin/activate"
    echo "/lb/project/mugqic/projects/MOH/moh_automation/moh_automation_main/transfer2json.py --input $TEMP/$LISTFILE --output /lb/project/mugqic/projects/MOH/Transfer_json --operation_cmd_line \"globus transfer --submission-id $sub_id --label $label --batch $TEMP/$LISTFILE $ABA_EP $BEL_EP\""
    echo "pt-cli ingest transfer --input-json "
else
    echo "$task_id failed!"
fi