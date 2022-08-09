#!/usr/bin/env python3
import glob
import sys
import os
from datetime import date
from  DB_OPS import create_connection
#,extract_sample_metrics,extract_sample_details,extract_fileloc_field
#update_metrics_db,extract_sample_details,extract_fileloc_details,extract_timestamp_details,update_timestamp_details,update_fileloc_details,extract_sample_names


def main():
    #Change these after Testing
    #file locations.
    Output_RR = "/home/dipop/MOH/TEST/MAIN/raw_reads"
    Work_dir = "/home/dipop/MOH/TEST/MAIN"
    Input_RR = "/home/dipop/MOH/TEST/raw_reads"
    SQL_DB = "/home/dipop/MOH/TEST/TEST.db"

    #Connect to the db
    connection = create_connection(SQL_DB)

    #Populate the lists
    RNA_Samples = []
    Tumour_Samples = []
    Normal_Samples = []
    Bad_Names = []
    



if __name__ == '__main__':
    main()
