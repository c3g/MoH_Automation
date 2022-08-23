# Marathon of Hope Automation #

This repository contains all of the scripts for the automation of the MOH projects
and integration of genpipes
## Abacus Scripts ##

### Activate_globus.sh ###
Activates the globus endpoints for Create_readset_transfer_BAMS.sh  
**Must be updated for individual users**  
Usage: `Activate_globus.sh`  

### Create_readset_transfer_BAMS.sh ###
Takes the full path of the read folder as an input and transfers the BAMs/fastqs with the 
MoH prefix.In addition it tranfers over the key run processing metrics and generates 
a log file forlater usage with the database  
**Globus fields must be updated for individual users**  
Usage: `Create_readset_transfer_BAMS.sh PATH_TO_RUN_FOLDER`  

## Beluga Scripts ##

### Metrics_Update.py ###
This script takes a single file and parses it to update the database with metrics from
run processing and from both DNA and RNA Genpipes. It compares the extracted values to
 known acceptable values and adds them to the KEY_METRICS table within the MoH database.  
Usage: `Metrics_Update.py` 

### MOH_Check_Progress.py ###
Parses the file structure of MAIN and looks for the output files produced from samples in
the Sample table. It then updates STATUS table with any progress. It queries for all files,
so any deliverables that are removed will result in an incomplete listing. In addition this 
script populates/updates the Timestamps and File_Locations tables.   
Usage: `MOH_Check_Progress.py`  

### DB_OPS.py ###
This contains all of the functions necessary for interacting with a file. Required for 
other scripts.  

### MOH_ln_output.py ###
Takes a single sample input and hard links all of the delivarable files within a structured
directory with their final names. In addition it constructs the custom readme and log files.
All data is taken from the database.   
Usage: `MOH_ln_output.py Sample`  
**Currently has a coverage cutoff implemented**  

### Create_CSVs.sh ###
Dumps the tables within the database as csv's within the CSV folder of DATABASE.   
Usage: `Create_CSVs.sh`  

### Generate_pairs_readset.py ###
Searches the temporary raw_reads folder for matching DNA pairs or RNA. It moves the
fastqs/BAMs and creates the readset and pairs files in preparation for Setup_run_XXX.sh. 
In addition, it populates the Samples table. It will not move process any files that are
outside the naming convention and it will exit if it finds files with the same name at the 
destination.    
Usage `Generate_pairs_readset.py DNA` or `Generate_pairs_readset.py RNA` 

