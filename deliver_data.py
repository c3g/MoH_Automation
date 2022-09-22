#!/usr/bin/env python3

#import glob
import os
import datetime
import pandas as pd
from  moh_resources import create_connection,extract_sample_metrics,extract_fileloc_field,extract_value,extract_sample_names, PatientData, getime

def main():
    connection = create_connection(r"/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db")
    #Test String
    #ALL_Samples = ["MoHQ-CM-3-1-17057-1D","MoHQ-CM-3-2-61285-1D","MoHQ-CM-1-9-2804-1-D"]
    samples_list = extract_sample_names(connection)
    for sample in samples_list:
        patient = PatientData(connection, sample)
        print(patient.patient)
        #check if samples reach the threashold for delivery
        #Check that dna_t dedup coverage is over 80
        #Check that dna_n dedup coverage is over 30
        #Check that processing is complete
        #Check that rna spots is over 100000000
        dna = False
        rna = False
        if  extract_sample_metrics(patient.conn,patient.dna_n,"WGS_Dedup_Coverage") or extract_value(patient.conn,"STATUS",patient.patient,"Tumour_Pair_Complete"):
            dna = False
        elif float(extract_sample_metrics(patient.conn,patient.dna_n,"WGS_Dedup_Coverage")) > 30 and float(extract_sample_metrics(patient.conn,patient.dna_t,"WGS_Dedup_Coverage")) >80:
            dna = True
        if  extract_sample_metrics(patient.conn,patient.rna,"WTS_Clusters") or extract_value(patient.conn,"STATUS",patient.patient,"rna_Complete"):
            rna = False
        elif float(extract_sample_metrics(patient.conn,patient.rna,"WTS_Clusters")) > 100000000:
            rna = True

        #Folders used for Delivery
        globus_share_folder = '/lustre03/project/rrg-bourqueg-ad/C3G/projects/GLOBUS_SHARE/MOH' # Output Folder
        globus_share_folder_patient = os.path.join(globus_share_folder, patient.institution, patient.cohort, patient.patient_corrected)
        #contains Warnings.txt Readme.txt Log.txt and all subfolders
        globus_share_raw_data_folder = os.path.join(globus_share_folder_patient, "raw_data")
        #contains raw bams and fastqs
        globus_share_variants_folder = os.path.join(globus_share_folder_patient, "variants")
        #contains all variants the subfolder
        globus_share_variants_callers_folder = os.path.join(globus_share_variants_folder, "caller_vcfs")
        #contains all the vcfs from the callers
        globus_share_alignment_folder = os.path.join(globus_share_folder_patient, "alignment")
        #contains the analysis bams
        globus_share_parameters_folder = os.path.join(globus_share_folder_patient, "paramaters")
        #contains the ini files.
        globus_share_tracks_folder = os.path.join(globus_share_folder_patient, "tracks")
        #contains the big wig tracks
        globus_share_reports_folder = os.path.join(globus_share_folder_patient, "reports")
        #contains the big wig tracks

        #I AM SO SORRY. I was getting frustrated and I got lazy and made two global varaiables. I have sinned.
        #UPDATED keeps track of if we need to update the metrics and warnings file
        #old_log stores the data in the log file to see if things need updating.
        # global UPDATED
        updated = False
        # global Old_log
        previous_log = {}
        log_file = os.path.join(globus_share_folder_patient + "log.txt")
        #See if the directory is created and if so check for file updates.
        if os.path.isdir(globus_share_folder_patient):
            #Load the log file into a dictionary for checking if updates are necessary.
            with open(log_file, "r", encoding="utf-8") as log_in:
                print ("logging")
                for line in log_in:
                    line.rstrip()
                    fields = line.split(",")
                    previous_log[fields[0]] = fields[1]

            # log = open(log_file, "a", encoding="utf-8")

        elif dna or rna:
            os.makedirs(globus_share_folder_patient)
            os.makedirs(globus_share_raw_data_folder)
            os.makedirs(globus_share_variants_folder)
            os.makedirs(globus_share_variants_callers_folder)
            os.makedirs(globus_share_alignment_folder)
            os.makedirs(globus_share_parameters_folder)
            os.makedirs(globus_share_tracks_folder)
            os.makedirs(globus_share_reports_folder)

            #Populate the general files
            with open(log_file, "w", encoding="utf-8") as log_in:
            # log = open(log_file, "w", encoding="utf-8")
                log_in.write("File,File_Creation_Date,Date_Added,Details\n")

            generate_readme(
                globus_share_folder_patient,
                patient.patient_corrected,
                patient.dna_n,
                patient.dna_t,
                patient.rna
                )
            log_new("Readme.txt",log_file,None,"Created")
        else:
            continue

        #Populate dna data
        if dna:
            updated = get_link_log(
                "dna_raw_bam_n",
                globus_share_raw_data_folder,
                ".bam",
                patient.dna_n,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_raw_bam_t",
                globus_share_raw_data_folder,
                ".bam",
                patient.dna_t,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

            updated = get_link_log(
                "dna_germline_vcf",
                globus_share_variants_folder,
                ".ensemble.germline.vt.annot.vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_somatic_vcf",
                globus_share_variants_folder,
                ".ensemble.somatic.vt.annot.vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_mutect2_germline_vcf",
                globus_share_variants_callers_folder,
                ".mutect2.germline.vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_mutect2_somatic_vcf",
                globus_share_variants_callers_folder,
                ".mutect2.somatic.vt.vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_strelka2_germline_vcf",
                globus_share_variants_callers_folder,
                ".strelka2.germline.vt.vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_strelka2_somatic_vcf",
                globus_share_variants_callers_folder,
                ".strelka2.somatic.vt.vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_vardict_germline_vcf",
                globus_share_variants_callers_folder,
                ".vardict.germline.vt.vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_vardict_somatic_vcf",
                globus_share_variants_callers_folder,
                ".vardict.somatic.vt.vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_varscan2_germline_vcf",
                globus_share_variants_callers_folder,
                ".varscan2.germline.vt.vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_varscan2_somatic_vcf",
                globus_share_variants_callers_folder,
                ".varscan2.somatic.vt.vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

            updated = get_link_log(
                "dna_final_bam_n",
                globus_share_alignment_folder,
                ".bam",
                patient.dna_n,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_final_bam_t",
                globus_share_alignment_folder,
                ".bam",
                patient.dna_t,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

            updated = get_link_log(
                "dna_multiqc_report",
                globus_share_reports_folder,
                "_D.multiqc.html",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "dna_pcgr_report",
                globus_share_reports_folder,
                ".pcgr.html",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

            updated = get_link_log(
                "dna_tumour_pair_ini",
                globus_share_parameters_folder,
                ".tp.ini",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

        if rna:
            print ("rna_TEST")
            updated = get_link_log(
                "rna_raw_fastq1",
                globus_share_raw_data_folder,
                "_R1.fastq.gz",
                patient.rna,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "rna_raw_fastq2",
                globus_share_raw_data_folder,
                "_R2.fastq.gz",
                patient.rna,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

            updated = get_link_log(
                "rna_vcf",
                globus_share_variants_folder,
                ".rna.hc.vcf.gz",
                patient.rna,connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

            updated = get_link_log(
                "rna_final_bam_expression",
                globus_share_alignment_folder,
                ".expression.cram",
                patient.rna,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "rna_final_bam_variants",
                globus_share_alignment_folder,
                "RT.variants.bam",
                patient.rna,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

            updated = get_link_log(
                "rna_multiqc_report",
                globus_share_reports_folder,
                "_R.multiqc.html",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "rna_annofuse_tsv",
                globus_share_reports_folder,
                ".anno_fuse",
                patient.rna,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "rna_gridss_report",
                globus_share_reports_folder,
                ".gridss",
                patient.rna,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

            updated = get_link_log(
                "rna_abundance_ini",
                globus_share_parameters_folder,
                ".rna.expression.ini",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "rna_variants_ini",
                globus_share_parameters_folder,
                ".rna.variants.ini",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

            updated = get_link_log(
                "rna_forward_bigwig",
                globus_share_tracks_folder,
                ".forward.bw",
                patient.rna,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )
            updated = get_link_log(
                "rna_reverse_bigwig",
                globus_share_tracks_folder,
                ".reverse.bw",
                patient.rna,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

        if rna and dna:
            updated = get_link_log(
                "dna_final_vcf",
                globus_share_variants_folder,
                ".vcf.gz",
                patient.patient_corrected,
                connection,
                log_file,
                patient.patient,
                updated,
                previous_log
                )

        #If any updates were made, Delete the old warning file and populate a new one.
        if updated:
            if os.path.exists(globus_share_folder_patient + "Warnings.txt"):
                os.remove(globus_share_folder_patient + "Warnings.txt")
                os.remove(globus_share_folder_patient + "Warnings.txt")
                log_new("Warnings.txt", log_file, None," Updated")
                log_new("Key_metrics.csv", log_file, None, "Updated")
            else:
                log_new("Warnings.txt", log_file, None, "Created")
                log_new("Key_metrics.csv", log_file, None, "Updated")

            #add key metrics table for samples
            metrics = pd.read_sql_query(f'select * from KEY_METRICS where patient="{patient.dna_n}" or patient="{patient.dna_t}" or patient="{patient.rna}"', connection)
            metrics.to_csv(globus_share_reports_folder + patient.patient_corrected + ".Key_metrics.csv", index=False)

            #Add warnings file
            warnings = pd.read_sql_query(f'select patient,YELLOW_Flags,RED_Flags from KEY_METRICS where patient="{patient.dna_n}" or patient="{patient.dna_t}" or patient="{patient.rna}"', connection)
            warnings.to_csv(globus_share_folder_patient + "Warnings.txt", index=False)
            with open(globus_share_folder_patient + 'Warnings.txt', 'r+', encoding="utf-8") as file:
                content = file.read()
                file.seek(0)
                file.write("Below are three collumns, Yellow flags indicate values that may be troublesome while red flags indicate a point of failure. Data may be useable with these flags, but any red flaged data should be carefully considered. If nothing is present, this data exceeded all standards.\n Data will not be delivered for Red Flaged coverage at this time.\n" + content)
        # log.close()

#Key function. Basiclly updates the logs and links the files based on the database.
def get_link_log(column, location, suffix, corrected_name, connection, log, name, update_status, previous_log):
    data = extract_fileloc_field(connection, name, column)
    if data:
        new_time = getime(data)
        #Portion for updating files if they have been modified
        if os.path.exists(location + corrected_name + suffix):
            old_time = previous_log[corrected_name + suffix]
            if old_time != new_time:
                os.remove(location + corrected_name + suffix)
                os.link(data,location + corrected_name + suffix)
                log_new(corrected_name + suffix, log,new_time,  "Updated")
                update_status = True
        else:
            os.link(data,location + corrected_name + suffix)
            log_new(corrected_name + suffix, log, new_time, "File Added")
            update_status = True
    return update_status
# def getime(PATH):
#     date = datetime.datetime.fromtimestamp(os.path.getmtime(PATH))
#     return date.strftime("%Y/%m/%d")

# class PatientData:
#     def __init__(self, connection, sample):
#         data = []
#         self.conn = connection
#         data = extract_patient_details(connection, sample)
#         if len(data) <10:
#             raise Exception(f'No database entry for {sample}')
#         self.sample = data[0]
#         self.Sample_True = data[1]
#         self.institution = data[2]
#         self.cohort  = data[3]
#         self.dna_n  = data[4]
#         self.dna_n_True = data[5]
#         self.dna_t = data[6]
#         self.dna_t_True = data[7]
#         self.rna = data[8]
#         self.rna_True = data[9]
#     def __str__(self):
#         return str(self.__class__) + ": " + str(self.__dict__)


def log_new(file, log_file, file_date, message):
##############Update this!!!
    if not file_date:
        file_date = datetime.date.today()
        file_date = file_date.strftime("%Y/%m/%d")
    now = datetime.date.today()
    date_formatted = now.strftime("%Y/%m/%d")
    print(file + "," + file_date + "," + date_formatted + "," + f"{message}\n")
    data = file + "," + file_date + "," + date_formatted + "," + f"{message}\n"
    with open(log_file, "a", encoding="utf-8") as log_in:
        log_in.write(data)



def generate_readme(location, patient, dna_n, dna_t, rna):
    text = f"""This directory contains the delivered data for {patient} processed by the Canadian Centre for Computational Genomics.
The data will be updated as it becomes available and as such many files may be missing from rna or dna upon initial creation of this directory
Should you have concerns, questions, or suggestions, please contact the analysis team at moh-q@computationalgenomics.ca
Within this directory you will find the results of the analysis for a single patient contained in 6 subdirectories and three files:
Readme.txt      This file
log.txt         Log file containing the dates of transfers and if files have been updated
Warnings.txt    Contains details of any warnings and whether they caused a failure of this analysis.
raw_data/               Contains all of the bam's/fastqs from the sequencer
    {dna_n}.bam         Raw dna reads for the Normal sample
    {dna_t}.bam         Raw dna reads for the Tumor sample
    {rna}.fastq       Raw rna reads for the Tumor sample
variants/                                   Contains the vcfs related to variant calls
    {patient}.ensemble.germline.vt.annot.vcf.gz     Germline Variants found in any of the callers
    {patient}.ensemble.somatic.vt.annot.vcf.gz      Somatic Variants found in any of the callers
    {patient}.rna.hc.vcf.gz                         Variants found using rna sample
*   {patient}.vcf.gz                                Contains the results of all callers for both dna and rna
    variants/caller_vcfs/                           Contains the vcfs produced from individual callers on the dna samples
        {patient}.mutect2.germline.vcf.gz           Germline results for mutect2
        {patient}.mutect2.somatic.vt.vcf.gz         Somatic results for mutect2
        {patient}.strelka2.germline.vt.vcf.gz       Germline results for strelka2
        {patient}.strelka2.somatic.vt.vcf.gz        Somatic results for strelka2
        {patient}.vardict.germline.vt.vcf.gz        Germline results for vardict
        {patient}.vardict.somatic.vt.vcf.gz         Somatic results for vardict
        {patient}.varscan2.germline.vt.vcf.gz       Germline results for varscan2
        {patient}.varscan2.somatic.vt.vcf.gz        Somatic results for varscan2
    
alignment/              Contains the alignment data for each sample
    {dna_n}.bam                 Alignment of normal against the reference
    {dna_t}.bam                 Alignment of tumor against the reference
    {rna}.expression.bam      Alignment of tumor rna against the reference used in expression analysis
    {rna}.variants.bam        Alignment of tumor rna against the reference used in variants analysis
reports/                    Contains the reports for the experiment
    {patient}_D.multiqc.html    QC report for the dna analysis
    {patient}_R.multiqc.html    QC report for the rna analysis
    {patient}.pcgr.html         Personal Cancer Genome Reporter report
    {rna}.anno_fuse         Report for fusions detected using rna
*   {rna}.gridss            Not yet available
    {patient}.Key_metrics.csv   Metrics used to determine whether the analysis was successful
parameters/             Contains the records of all the Parameters used in the experiment
    tumorPair.config.trace.ini              Parameters used in the tumor pair analysis
    rnaseq.expression.config.trace.ini      Parameters used in the rna expression analysis
*   rnaseq.variants.config.trace.ini        Parameters used in the rna variant analysis
tracks/                  Big Wig tracks for rna expression results
    {rna}.forward.bw        Forward big wig track
    {rna}.reverse.bw        Reverse big wig track
    """
    with open(os.path.join(location, "Readme.txt") , "w", encoding="utf-8") as filename:
        filename.write(text)


if __name__ == '__main__':
    main()
    #Update db with the objects
