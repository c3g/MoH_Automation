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

