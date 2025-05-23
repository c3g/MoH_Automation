#!/bin/bash

THIS_SCRIPT=$(basename "$0")

usage() {
    echo "script usage: $THIS_SCRIPT -h [-s specimen_name] [-l location] [-e experiment_nucleic_acid_type]"
    echo "Usage:"
    echo " -h                                 Display this help message."
    echo " -s <specimen_name>                 Specimen name to be delivered /!\ MANDATORY /!\."
    echo " -l <location>                      Location of the data to be delivered /!\ MANDATORY /!\."
    echo " -e <experiment_nucleic_acid_type>  Experiment nucleic acid type to be delivered (either: DNA or RNA) /!\ MANDATORY /!\."
    exit 1
}

while getopts 'hs:l:e:' OPTION; do
    case "$OPTION" in
        s)
            specimen_name="$OPTARG"
            ;;
        l)
            location="$OPTARG"
            ;;
        e)
            experiment_nucleic_acid_type="$OPTARG"
            case "$experiment_nucleic_acid_type" in
                DNA | RNA)
                    ;;
                *)
                    echo -e "ERROR: Invalid experiment nucleic acid type: '$experiment_nucleic_acid_type'.\n"
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

# Source cluster
if [[ $location == "abacus"* ]]; then
    # Abacus MOH
    SRC_MOH="/lb/project/mugqic/projects/MOH"
    # Abacus delivery folder location
    SRC_DELIVERY="$SRC_MOH/DELIVERY"
elif [[ $location == "beluga"* ]]; then
    # Beluga MOH
    SRC_MOH="/project/def-c3g/MOH"
    # Beluga delivery folder location
    SRC_DELIVERY="$SRC_MOH/DELIVERY"
elif [[ $location == "cardinal"* ]]; then
    # Cardinal MOH
    SRC_MOH="/project/def-c3g/MOH"
    # Cardinal delivery folder location
    SRC_DELIVERY="$SRC_MOH/DELIVERY"
else
    echo "ERROR: Unknown cluster. Exiting."
    exit 1
fi

if [ -z "${MUGQIC_INSTALL_HOME:-}" ]; then
  export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
fi

# mandatory arguments
if [ ! "$specimen_name" ] || [ ! "$location" ] || [ ! "$experiment_nucleic_acid_type" ]; then
  echo -e "ERROR: Missing mandatory arguments -s and/or -l and/or -e.\n"
  usage
fi

module avail 2>&1 | grep -m 1 -q "mugqic"; greprc=$?
if ! [[ $greprc -eq 0 ]]; then
  module use "$MUGQIC_INSTALL_HOME/modulefiles"
fi


timestamp=$(date "+%Y-%m-%dT%H.%M.%S")
delivery_json=${SRC_DELIVERY}/${specimen_name}_${experiment_nucleic_acid_type}_${timestamp}.json
listfile=${SRC_DELIVERY}/Delivery_${specimen_name}_${experiment_nucleic_acid_type}_${timestamp}.list
transfer_json=${SRC_DELIVERY}/${specimen_name}_${experiment_nucleic_acid_type}_${timestamp}_transfer.json

# shellcheck disable=SC1091
unset PYTHONPATH
source $SRC_MOH/project_tracking_cli/venv/bin/activate
# shellcheck disable=SC2086
pt-cli digest delivery --specimen_name $specimen_name --endpoint $location --experiment_nucleic_acid_type $experiment_nucleic_acid_type -o $delivery_json
# shellcheck disable=SC2086
$SRC_MOH/moh_automation/moh_automation_main/bucket_delivery.py -i $delivery_json -l $listfile
status=$?

if [ $status -ne 0 ]; then
    echo "Bucket delivery failed $SRC_MOH/moh_automation/moh_automation_main/bucket_delivery.py -i $delivery_json -l $listfile. Exiting..."
fi

$SRC_MOH/moh_automation/moh_automation_main/transfer2json.py --input $listfile --output $transfer_json --source $location --destination sd4h --operation_cmd_line "$SRC_MOH/moh_automation/moh_automation_main/bucket_delivery.py -i $delivery_json -l $listfile"
# shellcheck disable=SC2086
pt-cli ingest transfer --input-json $transfer_json
for i in $(cat $listfile | awk '{print $1}'); do
    rm $i
    location_id=$(pt-cli getid location --endpoint $location --file_name $(basename $i))
    # shellcheck disable=SC2086
    pt-cli --data "{\"modification\": [{\"table\": \"location\", \"id\": [\"$location_id\"]}]}" delete
done
