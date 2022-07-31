#!/usr/bin/env python3
import sqlite3
from sqlite3 import Error
import os
import re


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    return conn

def extract_sample_metrics(conn,Sample,column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM KEY_METRICS WHERE Sample='{Sample}'")
    result = cur.fetchone()
    if result:
        return result
    else:
        return "NA"
    cur.row_factory = None

def extract_sample_field(conn,Sample,column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM Samples WHERE Sample='{Sample}'")
    result = cur.fetchone()
    if result:
        return result
    else:
        return "NA"
    cur.row_factory = None

def extract_fileloc_field(conn,Sample,column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM File_Locations WHERE Sample='{Sample}'")
    result = cur.fetchone()
    if result:
        return result
    else:
        return "NA"
    cur.row_factory = None

def update_metrics_db(conn,Sample,WGS_Bases_Over_Q30 = 'NA',WGS_Min_Aligned_Reads_Delivered = 'NA',WGS_Raw_Coverage = 'NA',WGS_Dedup_Coverage = 'NA',Median_Insert_Size = 'NA',WGS_Duplication_Rate = 'NA',WGS_Contamination = 'NA',WTS_Clusters = 'NA',WTS_Unique_Reads='NA',WTS_Exonic_Rate='NA', WTS_rRNA_contamination = 'NA',Concordance = 'NA',Purity = 'NA',Yellow_Flags = 'NA',Red_Flags = 'NA'):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM KEY_METRICS WHERE Sample='{Sample}'""")
    result = cur.fetchone()
    if result:
        cur.execute(f"DELETE FROM KEY_METRICS WHERE Sample='{Sample}'")
        cur.execute(f"INSERT INTO KEY_METRICS (Sample,WGS_Bases_Over_Q30,WGS_Min_Aligned_Reads_Delivered,WGS_Raw_Coverage,WGS_Dedup_Coverage,Median_Insert_Size,WGS_Duplication_Rate,WGS_Contamination,WTS_Clusters,WTS_Exonic_Rate,WTS_Unique_Reads,WTS_rRNA_contamination,Concordance,Purity,Yellow_Flags,Red_Flags) VALUES ('{Sample}','{WGS_Bases_Over_Q30}','{WGS_Min_Aligned_Reads_Delivered}','{WGS_Raw_Coverage}','{WGS_Dedup_Coverage}','{Median_Insert_Size}','{WGS_Duplication_Rate}','{WGS_Contamination}','{WTS_Clusters}','{WTS_Exonic_Rate}','{WTS_Unique_Reads}','{WTS_rRNA_contamination}','{Concordance}','{Purity}','{Yellow_Flags}','{Red_Flags}')")
    else:
        cur.execute(f"INSERT INTO KEY_METRICS (Sample,WGS_Bases_Over_Q30,WGS_Min_Aligned_Reads_Delivered,WGS_Raw_Coverage,WGS_Dedup_Coverage,Median_Insert_Size,WGS_Duplication_Rate,WGS_Contamination,WTS_Clusters,WTS_Exonic_Rate,WTS_Unique_Reads,WTS_rRNA_contamination,Concordance,Purity,Yellow_Flags,Red_Flags) VALUES ('{Sample}','{WGS_Bases_Over_Q30}','{WGS_Min_Aligned_Reads_Delivered}','{WGS_Raw_Coverage}','{WGS_Dedup_Coverage}','{Median_Insert_Size}','{WGS_Duplication_Rate}','{WGS_Contamination}','{WTS_Clusters}','{WTS_Exonic_Rate}','{WTS_Unique_Reads}','{WTS_rRNA_contamination}','{Concordance}','{Purity}','{Yellow_Flags}','{Red_Flags}')")

        ############THIS FUNCTION NEEDS TO BE UPDATED FOR DATABASE CHANGES############
def update_sample_db(conn,sample,Institution,Cohort,DNA_N = 'ND',DNA_T = 'ND',RNA = 'ND'):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM Samples WHERE Sample='{sample}'""")
    result = cur.fetchone()

    if result:
        flag = False
        if result[3].startswith('MoH'):
            DNA_N = result[3]
            flag = True
        if result[4].startswith('MoH'):
            DNA_T = result[4]
            flag = True
        if result[5].startswith('MoH'):
            RNA = result[5]
            flag = True
        if flag == True:
            cur.execute(f"DELETE FROM Samples WHERE Sample='{sample}'")
            cur.execute(f"INSERT INTO Samples (Sample,Institution,Cohort,DNA_N,DNA_T,RNA) VALUES ('{sample}','{Institution}','{Cohort}','{DNA_N}','{DNA_T}','{RNA}')")
    else:
        cur.execute(f"INSERT INTO Samples (Sample,Institution,Cohort,DNA_N,DNA_T,RNA) VALUES ('{sample}','{Institution}','{Cohort}','{DNA_N}','{DNA_T}','{RNA}')")

#Grabs all fields from Sample table and outputs an array of the data
def extract_sample_details(conn,sample_true):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM Samples WHERE Sample='{sample_true}'""")
    result = cur.fetchone()
    if result:
        #return an aray that gets turned into an object for all the data from a sample.
        return result 
    else:
        return 'NA'

def extract_sample_names(conn):
    # Create cursor object
    cur = conn.cursor()
    cur.execute(f"""SELECT Sample FROM Samples""")
    cur.row_factory = lambda cursor, row: row[0]
    result = cur.fetchall()
    cur.row_factory = None
    if result:
        return result 
    else:
        return []

def extract_fileloc_details(conn,sample):
    cur = conn.cursor()
    cur.execute(f"""SELECT * FROM File_Locations WHERE Sample='{sample}'""")
    result = cur.fetchone()
    if result:
        return result 
    else:
        result = ['NA'] * 50
        return result

def extract_timestamp_details(conn,sample):
    cur = conn.cursor()
    cur.execute(f"""SELECT * FROM Timestamps WHERE Sample='{sample}'""")
    result = cur.fetchone()
    if result:
        #return an aray that gets turned into an object for all the data from a sample.
        return result 
    else:
        result = ['NA'] * 50
        return result

def update_timestamp_details(Sample):
    cur = Sample.conn.cursor()
    cur.execute(f"""SELECT * FROM Timestamps WHERE Sample='{Sample.Sample}'""")
    result = cur.fetchone()
    if result:
        cur.execute(f"DELETE FROM Timestamps WHERE Sample='{Sample.Sample}'")
        cur.execute(f"INSERT INTO Timestamps (Sample,Run_Proc_BAM_DNA_T,Run_Proc_BAM_DNA_N,Beluga_BAM_DNA_T,Beluga_BAM_DNA_N,DNA_VCF_G,DNA_VCF_S,Mutect2_Somatic_vcf,Mutect2_Germline_vcf,strelka2_Germline_vcf,strelka2_Somatic_vcf,vardict_Germline_vcf,vardict_Somatic_vcf,varscan2_Germline_vcf,varscan2_Somatic_vcf,Final_VCF,Final_DNA_BAM_T,Final_DNA_BAM_N,DNA_MultiQC,PCGR,TP_ini,Run_Proc_fastq_1_RNA,Run_Proc_fastq_2_RNA,Beluga_fastq_1_RNA,Beluga_fastq_2_RNA,RNA_VCF,Final_RNA_BAM_expression,Final_RNA_BAM_variants,RNA_MultiQC,AnnoFuse,GRIDSS,RNA_Abundance,big_wig_tracks_F,big_wig_tracks_R,RNA_Abundance_ini,RNA_Variants_ini) VALUES ('{Sample.Sample}','{Sample.TS_Run_Proc_BAM_DNA_T}','{Sample.TS_Run_Proc_BAM_DNA_N}','{Sample.TS_Beluga_BAM_DNA_T}','{Sample.TS_Beluga_BAM_DNA_N}','{Sample.TS_DNA_VCF_G}','{Sample.TS_DNA_VCF_S}','{Sample.TS_Mutect2_Somatic_vcf}','{Sample.TS_Mutect2_Germline_vcf}','{Sample.TS_strelka2_Germline_vcf}','{Sample.TS_strelka2_Somatic_vcf}','{Sample.TS_vardict_Germline_vcf}','{Sample.TS_vardict_Somatic_vcf}','{Sample.TS_varscan2_Germline_vcf}','{Sample.TS_varscan2_Somatic_vcf}','{Sample.TS_Final_VCF}','{Sample.TS_Final_DNA_BAM_T}','{Sample.TS_Final_DNA_BAM_N}','{Sample.TS_DNA_MultiQC}','{Sample.TS_PCGR}','{Sample.TS_TP_ini}','{Sample.TS_Run_Proc_fastq_1_RNA}','{Sample.TS_Run_Proc_fastq_2_RNA}','{Sample.TS_Beluga_fastq_1_RNA}','{Sample.TS_Beluga_fastq_2_RNA}','{Sample.TS_RNA_VCF}','{Sample.TS_Final_RNA_BAM_expression}','{Sample.TS_Final_RNA_BAM_variants}','{Sample.TS_RNA_MultiQC}','{Sample.TS_AnnoFuse}','{Sample.TS_GRIDSS}','{Sample.TS_RNA_Abundance}','{Sample.TS_big_wig_tracks_F}','{Sample.TS_big_wig_tracks_R}','{Sample.TS_RNA_Abundance_ini}','{Sample.TS_RNA_Variants_ini}')")
    else:
        cur.execute(f"INSERT INTO Timestamps (Sample,Run_Proc_BAM_DNA_T,Run_Proc_BAM_DNA_N,Beluga_BAM_DNA_T,Beluga_BAM_DNA_N,DNA_VCF_G,DNA_VCF_S,Mutect2_Somatic_vcf,Mutect2_Germline_vcf,strelka2_Germline_vcf,strelka2_Somatic_vcf,vardict_Germline_vcf,vardict_Somatic_vcf,varscan2_Germline_vcf,varscan2_Somatic_vcf,Final_VCF,Final_DNA_BAM_T,Final_DNA_BAM_N,DNA_MultiQC,PCGR,TP_ini,Run_Proc_fastq_1_RNA,Run_Proc_fastq_2_RNA,Beluga_fastq_1_RNA,Beluga_fastq_2_RNA,RNA_VCF,Final_RNA_BAM_expression,Final_RNA_BAM_variants,RNA_MultiQC,AnnoFuse,GRIDSS,RNA_Abundance,big_wig_tracks_F,big_wig_tracks_R,RNA_Abundance_ini,RNA_Variants_ini) VALUES ('{Sample.Sample}','{Sample.TS_Run_Proc_BAM_DNA_T}','{Sample.TS_Run_Proc_BAM_DNA_N}','{Sample.TS_Beluga_BAM_DNA_T}','{Sample.TS_Beluga_BAM_DNA_N}','{Sample.TS_DNA_VCF_G}','{Sample.TS_DNA_VCF_S}','{Sample.TS_Mutect2_Somatic_vcf}','{Sample.TS_Mutect2_Germline_vcf}','{Sample.TS_strelka2_Germline_vcf}','{Sample.TS_strelka2_Somatic_vcf}','{Sample.TS_vardict_Germline_vcf}','{Sample.TS_vardict_Somatic_vcf}','{Sample.TS_varscan2_Germline_vcf}','{Sample.TS_varscan2_Somatic_vcf}','{Sample.TS_Final_VCF}','{Sample.TS_Final_DNA_BAM_T}','{Sample.TS_Final_DNA_BAM_N}','{Sample.TS_DNA_MultiQC}','{Sample.TS_PCGR}','{Sample.TS_TP_ini}','{Sample.TS_Run_Proc_fastq_1_RNA}','{Sample.TS_Run_Proc_fastq_2_RNA}','{Sample.TS_Beluga_fastq_1_RNA}','{Sample.TS_Beluga_fastq_2_RNA}','{Sample.TS_RNA_VCF}','{Sample.TS_Final_RNA_BAM_expression}','{Sample.TS_Final_RNA_BAM_variants}','{Sample.TS_RNA_MultiQC}','{Sample.TS_AnnoFuse}','{Sample.TS_GRIDSS}','{Sample.TS_RNA_Abundance}','{Sample.TS_big_wig_tracks_F}','{Sample.TS_big_wig_tracks_R}','{Sample.TS_RNA_Abundance_ini}','{Sample.TS_RNA_Variants_ini}')")

def update_fileloc_details(Sample):
    cur = Sample.conn.cursor()
    cur.execute(f"""SELECT * FROM File_Locations WHERE Sample='{Sample.Sample}'""")
    result = cur.fetchone()
    if result:
        cur.execute(f"DELETE FROM File_Locations WHERE Sample='{Sample.Sample}'")
        cur.execute(f"INSERT INTO File_Locations (Sample,Run_Proc_BAM_DNA_T,Run_Proc_BAM_DNA_N,Beluga_BAM_DNA_T,Beluga_BAM_DNA_N,DNA_VCF_G,DNA_VCF_S,Mutect2_Somatic_vcf,Mutect2_Germline_vcf,strelka2_Germline_vcf,strelka2_Somatic_vcf,vardict_Germline_vcf,vardict_Somatic_vcf,varscan2_Germline_vcf,varscan2_Somatic_vcf,Final_VCF,Final_DNA_BAM_T,Final_DNA_BAM_N,DNA_MultiQC,PCGR,TP_ini,Run_Proc_fastq_1_RNA,Run_Proc_fastq_2_RNA,Beluga_fastq_1_RNA,Beluga_fastq_2_RNA,RNA_VCF,Final_RNA_BAM_expression,Final_RNA_BAM_variants,RNA_MultiQC,AnnoFuse,GRIDSS,RNA_Abundance,big_wig_tracks_F,big_wig_tracks_R,RNA_Abundance_ini,RNA_Variants_ini) VALUES ('{Sample.Sample}','{Sample.Run_Proc_BAM_DNA_T}','{Sample.Run_Proc_BAM_DNA_N}','{Sample.Beluga_BAM_DNA_T}','{Sample.Beluga_BAM_DNA_N}','{Sample.DNA_VCF_G}','{Sample.DNA_VCF_S}','{Sample.Mutect2_Somatic_vcf}','{Sample.Mutect2_Germline_vcf}','{Sample.strelka2_Germline_vcf}','{Sample.strelka2_Somatic_vcf}','{Sample.vardict_Germline_vcf}','{Sample.vardict_Somatic_vcf}','{Sample.varscan2_Germline_vcf}','{Sample.varscan2_Somatic_vcf}','{Sample.Final_VCF}','{Sample.Final_DNA_BAM_T}','{Sample.Final_DNA_BAM_N}','{Sample.DNA_MultiQC}','{Sample.PCGR}','{Sample.TP_ini}','{Sample.Run_Proc_fastq_1_RNA}','{Sample.Run_Proc_fastq_2_RNA}','{Sample.Beluga_fastq_1_RNA}','{Sample.Beluga_fastq_2_RNA}','{Sample.RNA_VCF}','{Sample.Final_RNA_BAM_expression}','{Sample.Final_RNA_BAM_variants}','{Sample.RNA_MultiQC}','{Sample.AnnoFuse}','{Sample.GRIDSS}','{Sample.RNA_Abundance}','{Sample.big_wig_tracks_F}','{Sample.big_wig_tracks_R}','{Sample.RNA_Abundance_ini}','{Sample.RNA_Variants_ini}')")
    else:
        cur.execute(f"INSERT INTO File_Locations (Sample,Run_Proc_BAM_DNA_T,Run_Proc_BAM_DNA_N,Beluga_BAM_DNA_T,Beluga_BAM_DNA_N,DNA_VCF_G,DNA_VCF_S,Mutect2_Somatic_vcf,Mutect2_Germline_vcf,strelka2_Germline_vcf,strelka2_Somatic_vcf,vardict_Germline_vcf,vardict_Somatic_vcf,varscan2_Germline_vcf,varscan2_Somatic_vcf,Final_VCF,Final_DNA_BAM_T,Final_DNA_BAM_N,DNA_MultiQC,PCGR,TP_ini,Run_Proc_fastq_1_RNA,Run_Proc_fastq_2_RNA,Beluga_fastq_1_RNA,Beluga_fastq_2_RNA,RNA_VCF,Final_RNA_BAM_expression,Final_RNA_BAM_variants,RNA_MultiQC,AnnoFuse,GRIDSS,RNA_Abundance,big_wig_tracks_F,big_wig_tracks_R,RNA_Abundance_ini,RNA_Variants_ini) VALUES ('{Sample.Sample}','{Sample.Run_Proc_BAM_DNA_T}','{Sample.Run_Proc_BAM_DNA_N}','{Sample.Beluga_BAM_DNA_T}','{Sample.Beluga_BAM_DNA_N}','{Sample.DNA_VCF_G}','{Sample.DNA_VCF_S}','{Sample.Mutect2_Somatic_vcf}','{Sample.Mutect2_Germline_vcf}','{Sample.strelka2_Germline_vcf}','{Sample.strelka2_Somatic_vcf}','{Sample.vardict_Germline_vcf}','{Sample.vardict_Somatic_vcf}','{Sample.varscan2_Germline_vcf}','{Sample.varscan2_Somatic_vcf}','{Sample.Final_VCF}','{Sample.Final_DNA_BAM_T}','{Sample.Final_DNA_BAM_N}','{Sample.DNA_MultiQC}','{Sample.PCGR}','{Sample.TP_ini}','{Sample.Run_Proc_fastq_1_RNA}','{Sample.Run_Proc_fastq_2_RNA}','{Sample.Beluga_fastq_1_RNA}','{Sample.Beluga_fastq_2_RNA}','{Sample.RNA_VCF}','{Sample.Final_RNA_BAM_expression}','{Sample.Final_RNA_BAM_variants}','{Sample.RNA_MultiQC}','{Sample.AnnoFuse}','{Sample.GRIDSS}','{Sample.RNA_Abundance}','{Sample.big_wig_tracks_F}','{Sample.big_wig_tracks_R}','{Sample.RNA_Abundance_ini}','{Sample.RNA_Variants_ini}')")


def update_status_db(conn,Sample,DNA_N_Transfered,DNA_T_Transfered,Alignment,Variants,Reports,Tumour_Pair_Complete,RNA_Transfrered,RNA_Alignment_expression,RNA_Alignment_Variant,RNA_Reports,RNA_Complete):
    cur = conn.cursor()
    cur.execute(f"SELECT * FROM STATUS WHERE Sample='{Sample}'")
    result = cur.fetchone()
    if result:
        cur.execute(f"DELETE FROM STATUS WHERE Sample='{sample}'")
        cur.execute(f"INSERT INTO STATUS (Sample,DNA_N_Transfered,DNA_T_Transfered,Alignment,Variants,Reports,Tumour_Pair_Complete,RNA_Transfrered,RNA_Alignment_expression,RNA_Alignment_Variant,RNA_Reports,RNA_Complete) VALUES ('{Sample}','{DNA_N_Transfered}','{DNA_T_Transfered}','{Alignment}','{Variants}','{Reports}','{Tumour_Pair_Complete}','{RNA_Transfrered}','{RNA_Alignment_expression}','{RNA_Alignment_Variant}','{RNA_Reports}','{RNA_Complete}')")
    else:
        cur.execute(f"INSERT INTO STATUS (Sample,DNA_N_Transfered,DNA_T_Transfered,Alignment,Variants,Reports,Tumour_Pair_Complete,RNA_Transfrered,RNA_Alignment_expression,RNA_Alignment_Variant,RNA_Reports,RNA_Complete) VALUES ('{Sample}','{DNA_N_Transfered}','{DNA_T_Transfered}','{Alignment}','{Variants}','{Reports}','{Tumour_Pair_Complete}','{RNA_Transfrered}','{RNA_Alignment_expression}','{RNA_Alignment_Variant}','{RNA_Reports}','{RNA_Complete}')")






