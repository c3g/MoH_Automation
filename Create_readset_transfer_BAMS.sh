#!/usr/bin/env bash

#Temporary File Location
TEMP='/lb/project/mugqic/projects/MOH/TEMP'
TIMESTAMP=`date +%FT%H.%M.%S`
LOGFILE=$TIMESTAMP"_transfer_log.txt";
LISTFILE=$TIMESTAMP"_transfer_list.txt";
touch $TEMP"/"$LOGFILE;
touch $TEMP"/"$LISTFILE;
echo "Log file of transfer from abacus to Beluga">"$TEMP/"$LOGFILE;
#CWD=$(pwd)
CWD="/lb/project/mugqic/projects/MOH/"
#Location of processing data
Location="$1"
if [[ $Location == */ ]]
then 
	echo $Location
else
	Location=$Location"/"
fi
echo $Location
echo "Transfered From"$Location>>"$TEMP/"$LOGFILE;

#Location on Beluga. CURRENTLY VERY IMPORTANT. DO NOT CHANGE OR IT WILL BREAK THE DATABASE
#SERIOUSLY DON'T CHANGE IT.
#PLEASE DONT.
BEL_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/raw_reads/"
#Beluga log file location
BEL_LOG_LOC="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/"
#Beluga Endpoint
BEL_EP='278b9bfe-24da-11e9-9fa2-0a06afd4a22e'
#Narval Endpoint #not used
NAR_EP='a1713da6-098f-40e6-b3aa-034efe8b6e5b'
#abacus Endpoint
ABA_EP='6c66d53d-a79d-11e8-96fa-0a6d4e044368'

#Loop over BAMS
#This is specific for MoHQ named Bams!
if ls "$Location"Aligned*/*/*/*/MoHQ*.bam 1> /dev/null 2>&1; then
echo "Bams Transfered">>"$TEMP/"$LOGFILE;
for i in "$Location"Aligned*/*/*/*/MoHQ*.bam; do
	rootname=`echo "$i" | cut -d'/' -f11`
	#Make the oneliner readset file.
	touch "$TEMP/"$rootname"_readset.tsv";
	#Header could probably be removed, but currently in for testing purposes
	#And is it really that hard to not grab the top line?
	echo 'Sample	Readset	LibraryType	RunType	Run	Lane	Adapter1	Adapter2	QualityOffset	BED	FASTQ1	FASTQ2	BAM' > "$TEMP/"$rootname"_readset.tsv";
	DETAILS=$(cat "$Location"*run.csv|grep $rootname;)
	RUNTYPE=$( echo $DETAILS | awk -F "\"*,\"*" '{print $4}' )
	RUNID=$( echo $DETAILS | awk -F "\"*,\"*" '{print $2}' )
	LANE=$( echo $DETAILS | awk -F "\"*,\"*" '{print $3}' )
	ADAP1=$( echo $DETAILS | awk -F "\"*,\"*" '{print $28}' )
	ADAP2=$( echo $DETAILS | awk -F "\"*,\"*" '{print $29}' )
	QUAL_OF=33
	BED=""
	BAM="raw_reads/"$rootname"/"$rootname".bam"
	FASTQ1=""
	FASTQ2=""
	echo $rootname"	"$rootname"."$RUNID"_"$LANE"	"$RUNID"_"$LANE"	"$RUNTYPE"	"$RUNID"	"$LANE"	"$ADAP1"	"$ADAP2"	"$QUAL_OF"	"$BED"	"$FASTQ1"	"$FASTQ2"	"$BAM>>"$TEMP/"$rootname"_readset.tsv";
	#Transfer Readset File
##########ALTER THIS###########
#	globus transfer --notify off -s checksum $ABA_EP:$TEMP"/"$rootname"_readset.tsv" $BEL_EP:$BEL_LOC$rootname"/"$rootname"_readset.tsv";
	echo $TEMP"/"$rootname"_readset.tsv" $BEL_LOC$rootname"/"$rootname"_readset.tsv">>$TEMP"/"$LISTFILE;
##########ALTER THIS###########
	#Start the transfer
	#globus transfer --notify off -s size $ABA_EP:$i $BEL_EP:$BEL_LOC$rootname"/"$rootname".bam"
	echo $i $BEL_LOC$rootname"/"$rootname".bam">>$TEMP"/"$LISTFILE;
	echo "$i"","$rootname"/"$rootname".bam">>"$TEMP/"$LOGFILE;
	
done;
else
	echo "No BAM files found";
fi;

#RNA now or possibly both
if ls "$Location"Unaligned*/*/*/MoHQ*_R1_001.fastq.gz 1> /dev/null 2>&1; then
echo "fastq's Transfered">>"$TEMP/"$LOGFILE;
for i in "$Location"Unaligned*/*/*/MoHQ*_R1_001.fastq.gz; do
	j="${i/_R1_/_R2_}"
	rootname=`echo "$i" | cut -d'/' -f11`
	rootname="${rootname%_*}"	
	rootname="${rootname/Sample_/}"	
	echo $rootname
	#Make the oneliner readset file.
	touch "$TEMP/"$rootname"_readset.tsv";
	#Header could probably be removed, but currently in for testing purposes
	echo 'Sample	Readset	LibraryType	RunType	Run	Lane	Adapter1	Adapter2	QualityOffset	BED	FASTQ1	FASTQ2	BAM' > "$TEMP/"$rootname"_readset.tsv";
	DETAILS=$(cat "$Location"*run.csv|grep $rootname)
	RUNTYPE=$( echo $DETAILS | awk -F "\"*,\"*" '{print $4}' )
	RUNID=$( echo $DETAILS | awk -F "\"*,\"*" '{print $2}' )
	LANE=$( echo $DETAILS | awk -F "\"*,\"*" '{print $3}' )
	ADAP1=$( echo $DETAILS | awk -F "\"*,\"*" '{print $28}' )
	ADAP2=$( echo $DETAILS | awk -F "\"*,\"*" '{print $29}' )
	QUAL_OF=33
	BED=""
	BAM=""
	FASTQ1="raw_reads/"$rootname"/"$rootname"_R1.fastq.gz"
	FASTQ2="raw_reads/"$rootname"/"$rootname"_R2.fastq.gz"
	echo $rootname"	"$rootname"."$RUNID"_"$LANE"	"$RUNID"_"$LANE"	"$RUNTYPE"	"$RUNID"	"$LANE"	"$ADAP1"	"$ADAP2"	"$QUAL_OF"	"$BED"	"$FASTQ1"	"$FASTQ2"	"$BAM>>"$TEMP/"$rootname"_readset.tsv";

##########ALTER THIS###########
	#Transfer Readset File
#	globus transfer --notify off -s checksum $ABA_EP:$TEMP"/"$rootname"_readset.tsv" $BEL_EP:$BEL_LOC$rootname"/"$rootname"_readset.tsv";
	echo $TEMP"/"$rootname"_readset.tsv" $BEL_LOC$rootname"/"$rootname"_readset.tsv">>$TEMP"/"$LISTFILE;
##########ALTER THIS###########
	#Start the transfer
	#globus transfer --notify off -s size $ABA_EP:$i $BEL_EP:$BEL_LOC$rootname"/"$rootname"_R1.fastq.gz"
	#globus transfer --notify off -s size $ABA_EP:$j $BEL_EP:$BEL_LOC$rootname"/"$rootname"_R2.fastq.gz"
	echo $i $BEL_LOC$rootname"/"$rootname"_R1.fastq.gz">>$TEMP"/"$LISTFILE;
	echo $j $BEL_LOC$rootname"/"$rootname"_R2.fastq.gz">>$TEMP"/"$LISTFILE;

	echo "$i"","$rootname"/"$rootname"_R1.fastq.gz">>"TEMP/"$LOGFILE;
	echo "$j"","$rootname"/"$rootname"_R2.fastq.gz">>"TEMP/"$LOGFILE;
done;
else
	echo "No fastq files found";
fi;
##########ALTER THIS###########
#globus transfer --notify off -s size $ABA_EP:$CWD"TEMP/"$LOGFILE $BEL_EP:$BEL_LOG_LOC$LOGFILE
echo $TEMP"/"$LOGFILE $BEL_LOG_LOC$LOGFILE>>$TEMP"/"$LISTFILE;

#transfer over the metrics
MET_LOC=$(ls $Location/*-novaseq-run.align_bwa_mem.csv)
F_NAME=${MET_LOC##*/}
##########ALTER THIS###########
#globus transfer --notify off -s size $ABA_EP:$MET_LOC $BEL_EP:/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/$F_NAME
echo $MET_LOC '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/'$F_NAME
>>$TEMP"/"$LISTFILE;

#Start the batch transfer.
globus transfer --batch $TEMP"/"$LISTFILE $ABA_EP $BEL_EP

