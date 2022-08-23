#!/usr/bin/env python3
#import glob
import sys
import os
import datetime
import pandas as pd
from  DB_OPS import create_connection,extract_sample_metrics,extract_sample_details,extract_fileloc_field,extract_value,extract_sample_names

def main():
    connection = create_connection(r"/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db")
    #Test String
    #ALL_Samples = ["MoHQ-CM-3-1-17057-1D","MoHQ-CM-3-2-61285-1D","MoHQ-CM-1-9-2804-1-D"] 
    ALL_Samples = extract_sample_names(connection)
    for smple in ALL_Samples: 
        Sample = SampleData(connection,smple)
        print (Sample.Sample)
        #check if samples reach the threashold for delivery
        #Check that DNA_T dedup coverage is over 80
        #Check that DNA_N dedup coverage is over 30
        #Check that processing is complete
        #Check that RNA spots is over 100000000
        DNA = False
        RNA = False
        if  extract_sample_metrics(Sample.conn,Sample.DNA_N,"WGS_Dedup_Coverage") == "NA" or extract_value(Sample.conn,"STATUS",Sample.Sample,"Tumour_Pair_Complete") == "NA":
            DNA = False
        elif float(extract_sample_metrics(Sample.conn,Sample.DNA_N,"WGS_Dedup_Coverage")) >30 and float(extract_sample_metrics(Sample.conn,Sample.DNA_T,"WGS_Dedup_Coverage")) >80:
            DNA = True
        if  extract_sample_metrics(Sample.conn,Sample.RNA,"WTS_Clusters") == "NA" or extract_value(Sample.conn,"STATUS",Sample.Sample,"RNA_Complete") == "NA":
            RNA = False 
        elif float(extract_sample_metrics(Sample.conn,Sample.RNA,"WTS_Clusters")) >100000000:
            RNA = True

        #Folders used for Delivery
        Base_Folder = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/'   #Base Folder
        Out_Folder = '/lustre03/project/rrg-bourqueg-ad/C3G/projects/GLOBUS_SHARE/MOH/' # Output Folder
        Out_Folder = Out_Folder + Sample.Institution + "/" + Sample.Cohort + "/" + Sample.Sample_True + "/"
        #contains Warnings.txt Readme.txt Log.txt and all subfolders
        Raw_Folder = Out_Folder + "raw_data/"
        #contains raw bams and fastqs
        Var_Folder = Out_Folder + "variants/"
        #contains all variants the subfolder 
        Cal_Folder = Var_Folder + "caller_vcfs/"
        #contains all the vcfs from the callers
        Align_Folder = Out_Folder + "alignment/"
        #contains the analysis bams
        Param_Folder = Out_Folder + "paramaters/"
        #contains the ini files.
        Tracks_Folder= Out_Folder + "tracks/"
        #contains the big wig tracks
        Reports_Folder= Out_Folder + "reports/"
        #contains the big wig tracks

        #I AM SO SORRY. I was getting frustrated and I got lazy and made two global varaiables. I have sinned. 
        #UPDATED keeps track of if we need to update the metrics and warnings file
        #old_log stores the data in the log file to see if things need updating.
        global UPDATED 
        UPDATED = False 
        global Old_log
        Old_log = {}
        #See if the directory is created and if so check for file updates.
        if os.path.isdir(Out_Folder):
            #Load the log file into a dictionary for checking if updates are necessary.
            with open(Out_Folder + "log.txt") as log_in:
                print ("logging")
                for line in log_in:
                    line.rstrip()
                    fields = line.split(",")
                    Old_log[fields[0]] = fields[1]
            log_in.close()
            log = open(Out_Folder + "log.txt", "a") 

        elif DNA == True or RNA == True:
            os.makedirs(Out_Folder)
            os.makedirs(Raw_Folder)
            os.makedirs(Var_Folder)
            os.makedirs(Cal_Folder)
            os.makedirs(Align_Folder)
            os.makedirs(Param_Folder)
            os.makedirs(Tracks_Folder)
            os.makedirs(Reports_Folder)

            #Populate the general files 
            log = open(Out_Folder + "log.txt", "w") 
            log.write("File,File_Creation_Date,Date_Added,Details\n")

            generate_readme(Out_Folder,Sample.Sample_True,Sample.DNA_N,Sample.DNA_T,Sample.RNA)
            log_new("Readme.txt",log,None,"Created")
        else:
            continue



#Populate DNA data
        if DNA == True:
            get_link_log("Beluga_BAM_DNA_N",Raw_Folder,".bam",Sample.DNA_N,connection,log,Sample.Sample)
            get_link_log("Beluga_BAM_DNA_T",Raw_Folder,".bam",Sample.DNA_T,connection,log,Sample.Sample)
            
            get_link_log("DNA_VCF_G",Var_Folder,".ensemble.germline.vt.annot.vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("DNA_VCF_S",Var_Folder,".ensemble.somatic.vt.annot.vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("Mutect2_Germline_vcf",Cal_Folder,".mutect2.germline.vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("Mutect2_Somatic_vcf",Cal_Folder,".mutect2.somatic.vt.vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("strelka2_Germline_vcf",Cal_Folder,".strelka2.germline.vt.vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("strelka2_Somatic_vcf",Cal_Folder,".strelka2.somatic.vt.vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("vardict_Germline_vcf",Cal_Folder,".vardict.germline.vt.vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("vardict_Somatic_vcf",Cal_Folder,".vardict.somatic.vt.vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("varscan2_Germline_vcf",Cal_Folder,".varscan2.germline.vt.vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("varscan2_Somatic_vcf",Cal_Folder,".varscan2.somatic.vt.vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)

            get_link_log("Final_DNA_BAM_N",Align_Folder,".bam",Sample.DNA_N,connection,log,Sample.Sample)
            get_link_log("Final_DNA_BAM_T",Align_Folder,".bam",Sample.DNA_T,connection,log,Sample.Sample)

            get_link_log("DNA_MultiQC",Reports_Folder,"_D.multiqc.html",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("PCGR",Reports_Folder,".pcgr.html",Sample.Sample_True,connection,log,Sample.Sample)

            get_link_log("TP_ini",Param_Folder,".tp.ini",Sample.Sample_True,connection,log,Sample.Sample)

        if RNA == True:
            print ("RNA_TEST")
            get_link_log("Beluga_fastq_1_RNA",Raw_Folder,"_R1.fastq.gz",Sample.RNA,connection,log,Sample.Sample)
            get_link_log("Beluga_fastq_2_RNA",Raw_Folder,"_R2.fastq.gz",Sample.RNA,connection,log,Sample.Sample)

            get_link_log("RNA_VCF",Var_Folder,".rna.hc.vcf.gz",Sample.RNA,connection,log,Sample.Sample)

            get_link_log("Final_RNA_BAM_expression",Align_Folder,".expression.cram",Sample.RNA,connection,log,Sample.Sample)
            get_link_log("Final_RNA_BAM_variants",Align_Folder,"RT.variants.bam",Sample.RNA,connection,log,Sample.Sample)

            get_link_log("RNA_MultiQC",Reports_Folder,"_R.multiqc.html",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("AnnoFuse",Reports_Folder,".anno_fuse",Sample.RNA,connection,log,Sample.Sample)
            get_link_log("GRIDSS",Reports_Folder,".gridss",Sample.RNA,connection,log,Sample.Sample)

            get_link_log("RNA_Abundance_ini",Param_Folder,".RNA.expression.ini",Sample.Sample_True,connection,log,Sample.Sample)
            get_link_log("RNA_Variants_ini",Param_Folder,".RNA.variants.ini",Sample.Sample_True,connection,log,Sample.Sample)

            get_link_log("big_wig_tracks_F",Tracks_Folder,".forward.bw",Sample.RNA,connection,log,Sample.Sample)
            get_link_log("big_wig_tracks_R",Tracks_Folder,".reverse.bw",Sample.RNA,connection,log,Sample.Sample)

        if RNA == True and DNA == True:
            get_link_log("Final_VCF",Var_Folder,".vcf.gz",Sample.Sample_True,connection,log,Sample.Sample)

        
        #If any updates were made, Delete the old warning file and populate a new one.
        if UPDATED == True:
            if os.path.exists(Out_Folder + "Warnings.txt"):
                os.remove(Out_Folder + "Warnings.txt")
                os.remove(Out_Folder + "Warnings.txt")
                log_new("Warnings.txt",log,None,"Updated")
                log_new("Key_metrics.csv",log,None,"Updated")
            else:
                log_new("Warnings.txt",log,None,"Created")
                log_new("Key_metrics.csv",log,None,"Updated")

            #add key metrics table for samples
            metrics = pd.read_sql_query(f'select * from KEY_METRICS where Sample="{Sample.DNA_N}" or Sample="{Sample.DNA_T}" or Sample="{Sample.RNA}"', connection) 
            metrics.to_csv(Reports_Folder + Sample.Sample_True + ".Key_metrics.csv", index=False)

            #Add warnings file
            Warnings = pd.read_sql_query(f'select Sample,YELLOW_Flags,RED_Flags from KEY_METRICS where Sample="{Sample.DNA_N}" or Sample="{Sample.DNA_T}" or Sample="{Sample.RNA}"', connection) 
            Warnings.to_csv(Out_Folder + "Warnings.txt", index=False)
            with open(Out_Folder + 'Warnings.txt', 'r+') as file:
                content = file.read()
                file.seek(0)
                file.write("Below are three collumns, Yellow flags indicate values that may be troublesome while red flags indicate a point of failure. Data may be useable with these flags, but any red flaged data should be carefully considered. If nothing is present, this data exceeded all standards.\n Data will not be delivered for Red Flaged coverage at this time.\n" + content)
        log.close()

#Key function. Basiclly updates the logs and links the files based on the database.
def get_link_log(Column,location,suffix,True_Name,connection,log,Name):
    data = extract_fileloc_field(connection,Name,Column)
    if data != "NA":
        new_time = getime(data)
        #Portion for updating files if they have been modified
        if os.path.exists(location + True_Name + suffix):
            old_time = Old_log[True_Name + suffix] 
            if old_time != new_time:
                os.remove(location + True_Name + suffix)
                os.link(data,location + True_Name + suffix)
                log_new(True_Name + suffix, log,new_time,  "Updated")
                UPDATED = True
        else:
            os.link(data,location + True_Name + suffix)
            log_new(True_Name + suffix, log, new_time, "File Added")
            UPDATED = True

def getime(PATH):
    date = datetime.datetime.fromtimestamp(os.path.getmtime(PATH))
    return date.strftime("%Y/%m/%d")

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


def log_new(File,log,file_date,message):
##############Update this!!!
    if file_date == None:
        file_date = datetime.date.today()
        file_date = file_date.strftime("%Y/%m/%d")
    Now = datetime.date.today()
    Date = Now.strftime("%Y/%m/%d")
    print (File + "," + file_date + "," + Date + "," + f"{message}\n")
    data = File + "," + file_date + "," + Date + "," + f"{message}\n"
    log.write(data)



def generate_readme(Location,PATIENT,DNA_N,DNA_T,RNA):
    data = f"""This directory contains the delivered data for {PATIENT} processed by the Canadian Centre for Computational Genomics. 
The data will be updated as it becomes available and as such many files may be missing from RNA or DNA upon initial creation of this directory
Should you have concerns, questions, or suggestions, please contact the analysis team at moh-q@computationalgenomics.ca
Within this directory you will find the results of the analysis for a single patient contained in 6 subdirectories and three files:
Readme.txt      This file
log.txt         Log file containing the dates of transfers and if files have been updated
Warnings.txt    Contains details of any warnings and whether they caused a failure of this analysis.
raw_data/               Contains all of the bam's/fastqs from the sequencer
    {DNA_N}.bam         Raw DNA reads for the Normal sample
    {DNA_T}.bam         Raw DNA reads for the Tumor sample
    {RNA}.fastq       Raw RNA reads for the Tumor sample
variants/                                   Contains the vcfs related to variant calls
    {PATIENT}.ensemble.germline.vt.annot.vcf.gz     Germline Variants found in any of the callers
    {PATIENT}.ensemble.somatic.vt.annot.vcf.gz      Somatic Variants found in any of the callers
    {PATIENT}.rna.hc.vcf.gz                         Variants found using RNA sample
*   {PATIENT}.vcf.gz                                Contains the results of all callers for both DNA and RNA
    variants/caller_vcfs/                           Contains the vcfs produced from individual callers on the DNA samples
        {PATIENT}.mutect2.germline.vcf.gz           Germline results for mutect2
        {PATIENT}.mutect2.somatic.vt.vcf.gz         Somatic results for mutect2
        {PATIENT}.strelka2.germline.vt.vcf.gz       Germline results for strelka2
        {PATIENT}.strelka2.somatic.vt.vcf.gz        Somatic results for strelka2
        {PATIENT}.vardict.germline.vt.vcf.gz        Germline results for vardict
        {PATIENT}.vardict.somatic.vt.vcf.gz         Somatic results for vardict
        {PATIENT}.varscan2.germline.vt.vcf.gz       Germline results for varscan2
        {PATIENT}.varscan2.somatic.vt.vcf.gz        Somatic results for varscan2
    
alignment/              Contains the alignment data for each sample
    {DNA_N}.bam                 Alignment of normal against the reference
    {DNA_T}.bam                 Alignment of tumor against the reference
    {RNA}.expression.bam      Alignment of tumor RNA against the reference used in expression analysis
    {RNA}.variants.bam        Alignment of tumor RNA against the reference used in variants analysis
reports/                    Contains the reports for the experiment
    {PATIENT}_D.multiqc.html    QC report for the DNA analysis
    {PATIENT}_R.multiqc.html    QC report for the RNA analysis
    {PATIENT}.pcgr.html         Personal Cancer Genome Reporter report
    {RNA}.anno_fuse         Report for fusions detected using RNA
*   {RNA}.gridss            Not yet available
    {PATIENT}.Key_metrics.csv   Metrics used to determine whether the analysis was successful
parameters/             Contains the records of all the Parameters used in the experiment
    tumorPair.config.trace.ini              Parameters used in the tumor pair analysis
    RNAseq.expression.config.trace.ini      Parameters used in the RNA expression analysis
*   RNAseq.variants.config.trace.ini        Parameters used in the RNA variant analysis
tracks/                  Big Wig tracks for RNA expression results
    {RNA}.forward.bw        Forward big wig track
    {RNA}.reverse.bw        Reverse big wig track
    """
    f = open(Location + "Readme.txt" , "w") 
    f.write(data)
    f.close()



if __name__ == '__main__':
    main()
    #Update db with the objects
