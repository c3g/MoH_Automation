#!/bin/bash
set -eu -o pipefail

##################################################
# Setup module
module load mugqic/python/3.10.4

##################################################
# STEP 1: Update file paths and timestamps in DB
echo "Running MOH_Check_Progress.py ..."
./MOH_Check_Progress.py
echo -e "\n"

# Make sure previous job is done interacting with db
sleep 30

##################################################
# STEP 2: Update metrics in DB
echo "Running Metrics_Update.py ..."
./Metrics_Update.py
echo -e "\n"

# Make sure previous job is done interacting with db
sleep 30

##################################################
# STEP 3: Update CSV files after DB modification
echo "Running Create_CSVs.sh ..."
./Create_CSVs.sh
echo -e "Done.\n"

# Make sure previous job is done interacting with db
sleep 30

##################################################
# STEP 4: Deliver data
echo "Running MOH_ln_output.py ..."
./MOH_ln_output.py