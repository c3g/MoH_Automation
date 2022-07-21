#!/usr/bin/env python3
import glob
import sys
import re
from  DB_OPS import create_connection,extract_sample_metrics,extract_sample_details
#update_metrics_db,extract_sample_details,extract_fileloc_details,extract_timestamp_details,update_timestamp_details,update_fileloc_details,extract_sample_names
def main():
    Name = sys.argv[1]          #Current name of the sample as found in Samples table
    print(Name)
    Base_Folder = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/'   #Base Folder
    Out_Folder = '/lustre03/project/rrg-bourqueg-ad/C3G/projects/GLOBUS_SHARE/MOH/' # Output Folder
    connection = create_connection(r"/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db")
    Sample = SampleData(connection,Name)
    #check if samples reach the threashold for delivery
    #Check that DNA_T dedup coverage is over 80
    #Check that DNA_N dedup coverage is over 30
    #Check that RNA spots is over 100000000
    DNA = False
    RNA = False
    if  extract_sample_metrics(Sample.conn,Sample.DNA_N,"WGS_Dedup_Coverage") == "NA":
        DNA = False
    elif float(extract_sample_metrics(Sample.conn,Sample.DNA_N,"WGS_Dedup_Coverage")) >30 and float(extract_sample_metrics(Sample.conn,Sample.DNA_T,"WGS_Dedup_Coverage")) >80:
        DNA = True
    if  extract_sample_metrics(Sample.conn,Sample.DNA_N,"WTS_Clusters") == "NA":
        RNA = False 
    elif extract_sample_metrics("WTS_Clusters") >100000000:
        RNA = True
    print (DNA)
    print (RNA)

#See if the directory is created and if 




#Check if output directory exists, If not create it and a simple log.
#If it does check the file names individually before overwriting.
class SampleData:
    def __init__(self, connection, Sample):
        data = []
        self.conn = connection
        data = extract_sample_details(connection,Sample)
        if len(data) <10:
            raise Exception(f'No database entry for {Sample}')
        self.Sample = data[0]
        self.Sample_True = data[1]
        self.Institution = data[2]
        self.Cohort  = data[3]
        self.DNA_N  = data[4]
        self.DNA_N_True = data[5]
        self.DNA_T = data[6]
        self.DNA_T_True = data[7]
        self.RNA = data[8]
        self.RNA_True = data[9]
    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

if __name__ == '__main__':
    main()
    #Update db with the objects
