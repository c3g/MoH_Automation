# Marathon of Hope Automation #

This repository contains all of the scripts for the automation of the MOH projects
and integration of genpipes
## Abacus Scripts ##

### Activate_globus.sh ###
Activates the globus endpoints for Create_readset_transfer_BAMS.sh  
**Must be updated for individual users**  
Usage: `Activate_globus.sh`  

## Beluga Scripts ##


### Create_readset_transfer_BAMS.sh ###
Takes the full path of the read folder as an input and transfers the BAMs/fastqs with the 
MoH prefix.In addition it tranfers over the key run processing metrics and generates 
a log file forlater usage with the database  
**Globus fields must be updated for individual users**  
Usage: `Create_readset_transfer_BAMS.sh PATH_TO_RUN_FOLDER`  

### Metrics_Update.py ###
This script takes a single file and parses it to update the database with metrics from
run processing and from both DNA and RNA Genpipes. It compares the extracted values to
 known acceptable values and adds them to the KEY_METRICS table within the MoH database.  
**Update needed to extract data directly from Database rather than a file**  
Usage: `Metrics_Update.py` 

### MOH_Check_Progress.py ###
Parses Samples table for all available samples and updates the File_location & Timestamp
tables with the appropriate data.
Usage: `MOH_Check_Progress.py`  

### DB_OPS.py ###
This contains all of the functions necessary for interacting with a file. Required for 
other scripts.  

