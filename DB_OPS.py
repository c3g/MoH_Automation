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

# Replace the other functions with this one
def extract_value(conn, table, sample, column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM {table} WHERE Sample='{sample}'")
    result = cur.fetchone()
    if result:
        return result
    else:
        return "NA"
    cur.row_factory = None

def extract_patient_status(conn, patient, column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM status WHERE patient='{patient}'")
    result = cur.fetchone()
    if result:
        return result
    else:
        return "NA"
    cur.row_factory = None

def extract_sample_metrics(conn, sample, column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM KEY_METRICS WHERE Sample='{sample}'")
    result = cur.fetchone()
    if result:
        return result
    else:
        return "NA"
    cur.row_factory = None

def extract_sample_field(conn,sample,column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM Samples WHERE Sample='{sample}'")
    result = cur.fetchone()
    if result:
        return result
    else:
        return "NA"
    cur.row_factory = None

def extract_fileloc_field(conn, sample, column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM File_Locations WHERE Sample='{sample}'")
    result = cur.fetchone()
    if result:
        return result
    else:
        return "NA"
    cur.row_factory = None

def update_metrics_db(
    conn,
    sample,
    WGS_Bases_Over_Q30,
    WGS_Min_Aligned_Reads_Delivered,
    Raw_Mean_Coverage,
    WGS_Dedup_Coverage,
    Median_Insert_Size,
    Raw_Duplication_Rate,
    WGS_Contamination,
    Raw_Reads_Count,
    WTS_Aligned_Reads,
    WTS_Exonic_Rate,
    WTS_rRNA_contamination,
    Concordance,
    Purity,
    Flags,
    Fails,
    Raw_Median_Insert_Size,
    Raw_Mean_Insert_Size
    ):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM KEY_METRICS WHERE Sample='{sample}'""")
    result = cur.fetchone()
    if result:
        cur.execute(f"DELETE FROM KEY_METRICS WHERE Sample='{sample}'")
        cur.execute(f"INSERT INTO KEY_METRICS (Sample,WGS_Bases_Over_Q30,WGS_Min_Aligned_Reads_Delivered,Raw_Mean_Coverage,WGS_Dedup_Coverage,Median_Insert_Size,Raw_Duplication_Rate,WGS_Contamination,Raw_Reads_Count,WTS_Exonic_Rate,WTS_Aligned_Reads,WTS_rRNA_contamination,Concordance,Purity,Flags,Fails,Raw_Median_Insert_Size,Raw_Mean_Insert_Size) VALUES ('{sample}','{WGS_Bases_Over_Q30}','{WGS_Min_Aligned_Reads_Delivered}','{Raw_Mean_Coverage}','{WGS_Dedup_Coverage}','{Median_Insert_Size}','{Raw_Duplication_Rate}','{WGS_Contamination}','{Raw_Reads_Count}','{WTS_Exonic_Rate}','{WTS_Aligned_Reads}','{WTS_rRNA_contamination}','{Concordance}','{Purity}','{Flags}','{Fails}','{Raw_Median_Insert_Size}','{Raw_Mean_Insert_Size}')")
    else:
        cur.execute(f"INSERT INTO KEY_METRICS (Sample,WGS_Bases_Over_Q30,WGS_Min_Aligned_Reads_Delivered,Raw_Mean_Coverage,WGS_Dedup_Coverage,Median_Insert_Size,Raw_Duplication_Rate,WGS_Contamination,Raw_Reads_Count,WTS_Exonic_Rate,WTS_Aligned_Reads,WTS_rRNA_contamination,Concordance,Purity,Flags,Fails,Raw_Median_Insert_Size,Raw_Mean_Insert_Size) VALUES ('{sample}','{WGS_Bases_Over_Q30}','{WGS_Min_Aligned_Reads_Delivered}','{Raw_Mean_Coverage}','{WGS_Dedup_Coverage}','{Median_Insert_Size}','{Raw_Duplication_Rate}','{WGS_Contamination}','{Raw_Reads_Count}','{WTS_Exonic_Rate}','{WTS_Aligned_Reads}','{WTS_rRNA_contamination}','{Concordance}','{Purity}','{Flags}','{Fails}','{Raw_Median_Insert_Size}','{Raw_Mean_Insert_Size}')")

def Update_Samples_Table(conn, sample, sample_True, Instituion, Cohort, DNA_N, DNA_N_True, DNA_T, DNA_T_True, RNA, RNA_True):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM Samples WHERE Sample='{sample}'""")
    result = cur.fetchone()
    if result:
        if result[4].startswith('MoH') and DNA_N == 'NA':
            DNA_N = result[4]
            DNA_N_True = result[4]
        if result[6].startswith('MoH') and DNA_T == 'NA':
            DNA_T = result[6]
            DNA_T_True = result[6]
        if result[8].startswith('MoH') and RNA == 'NA':
            RNA = result[8]
            RNA_True = result[8]
        cur.execute(f"DELETE FROM Samples WHERE Sample='{sample}'")
        cur.execute(f"INSERT INTO Samples (Sample,Sample_True,Instituion,Cohort,DNA_N,DNA_N_True,DNA_T,DNA_T_True,RNA,RNA_True) VALUES ('{sample}','{sample}','{Instituion}','{Cohort}','{DNA_N}','{DNA_N_True}','{DNA_T}','{DNA_T_True}','{RNA}','{RNA_True}')")
    else:
        cur.execute(f"INSERT INTO Samples (Sample,Sample_True,Instituion,Cohort,DNA_N,DNA_N_True,DNA_T,DNA_T_True,RNA,RNA_True) VALUES ('{sample}','{sample}','{Instituion}','{Cohort}','{DNA_N}','{DNA_N_True}','{DNA_T}','{DNA_T_True}','{RNA}','{RNA_True}')")

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

def update_timestamp_details(sample):
    cur = sample.conn.cursor()
    cur.execute(f"""SELECT * FROM Timestamps WHERE Sample='{sample.sample}'""")
    result = cur.fetchone()
    if result:
        cur.execute(f"DELETE FROM Timestamps WHERE Sample='{sample.sample}'")
        cur.execute(f"INSERT INTO Timestamps (Sample,Run_Proc_BAM_DNA_T,Run_Proc_BAM_DNA_N,Beluga_BAM_DNA_T,Beluga_BAM_DNA_N,DNA_VCF_G,DNA_VCF_S,Mutect2_Somatic_vcf,Mutect2_Germline_vcf,strelka2_Germline_vcf,strelka2_Somatic_vcf,vardict_Germline_vcf,vardict_Somatic_vcf,varscan2_Germline_vcf,varscan2_Somatic_vcf,cnvkit_vcf,Final_VCF,Final_DNA_BAM_T,Final_DNA_BAM_N,DNA_MultiQC,pcgr_report,pcgr_maf,pcgr_snvs_indels,pcgr_cna_segments,TP_ini,Run_Proc_fastq_1_RNA,Run_Proc_fastq_2_RNA,Beluga_fastq_1_RNA,Beluga_fastq_2_RNA,RNA_VCF,Final_RNA_BAM_variants,RNA_MultiQC,AnnoFuse,GRIDSS,RNA_Abundance,big_wig_tracks_F,big_wig_tracks_R,RNA_Abundance_ini,RNA_Variants_ini) VALUES ('{sample.sample}','{sample.ts_run_proc_bam_dna_t}','{sample.ts_run_proc_bam_dna_n}','{sample.ts_beluga_bam_dna_t}','{sample.ts_beluga_bam_dna_n}','{sample.ts_dna_vcf_g}','{sample.ts_dna_vcf_s}','{sample.ts_mutect2_somatic_vcf}','{sample.ts_mutect2_germline_vcf}','{sample.ts_strelka2_germline_vcf}','{sample.ts_strelka2_somatic_vcf}','{sample.ts_vardict_germline_vcf}','{sample.ts_vardict_somatic_vcf}','{sample.ts_varscan2_germline_vcf}','{sample.ts_varscan2_somatic_vcf}','{sample.ts_cnvkit_vcf}','{sample.ts_final_vcf}','{sample.ts_final_dna_bam_t}','{sample.ts_final_dna_bam_n}','{sample.ts_dna_multiqc}','{sample.ts_pcgr_report}','{sample.ts_pcgr_maf}','{sample.ts_pcgr_snvs_indels}','{sample.ts_pcgr_cna_segments}','{sample.ts_tp_ini}','{sample.ts_run_proc_fastq_1_rna}','{sample.ts_run_proc_fastq_2_rna}','{sample.ts_beluga_fastq_1_rna}','{sample.ts_beluga_fastq_2_rna}','{sample.ts_rna_vcf}','{sample.ts_final_rna_bam_variants}','{sample.ts_rna_multiqc}','{sample.ts_annofuse}','{sample.ts_gridss}','{sample.ts_rna_abundance}','{sample.ts_big_wig_tracks_f}','{sample.ts_big_wig_tracks_r}','{sample.ts_rna_abundance_ini}','{sample.ts_rna_variants_ini}')")
    else:
        cur.execute(f"INSERT INTO Timestamps (Sample,Run_Proc_BAM_DNA_T,Run_Proc_BAM_DNA_N,Beluga_BAM_DNA_T,Beluga_BAM_DNA_N,DNA_VCF_G,DNA_VCF_S,Mutect2_Somatic_vcf,Mutect2_Germline_vcf,strelka2_Germline_vcf,strelka2_Somatic_vcf,vardict_Germline_vcf,vardict_Somatic_vcf,varscan2_Germline_vcf,varscan2_Somatic_vcf,cnvkit_vcf,Final_VCF,Final_DNA_BAM_T,Final_DNA_BAM_N,DNA_MultiQC,pcgr_report,pcgr_maf,pcgr_snvs_indels,pcgr_cna_segments,TP_ini,Run_Proc_fastq_1_RNA,Run_Proc_fastq_2_RNA,Beluga_fastq_1_RNA,Beluga_fastq_2_RNA,RNA_VCF,Final_RNA_BAM_variants,RNA_MultiQC,AnnoFuse,GRIDSS,RNA_Abundance,big_wig_tracks_F,big_wig_tracks_R,RNA_Abundance_ini,RNA_Variants_ini) VALUES ('{sample.sample}','{sample.ts_run_proc_bam_dna_t}','{sample.ts_run_proc_bam_dna_n}','{sample.ts_beluga_bam_dna_t}','{sample.ts_beluga_bam_dna_n}','{sample.ts_dna_vcf_g}','{sample.ts_dna_vcf_s}','{sample.ts_mutect2_somatic_vcf}','{sample.ts_mutect2_germline_vcf}','{sample.ts_strelka2_germline_vcf}','{sample.ts_strelka2_somatic_vcf}','{sample.ts_vardict_germline_vcf}','{sample.ts_vardict_somatic_vcf}','{sample.ts_varscan2_germline_vcf}','{sample.ts_varscan2_somatic_vcf}','{sample.ts_cnvkit_vcf}','{sample.ts_final_vcf}','{sample.ts_final_dna_bam_t}','{sample.ts_final_dna_bam_n}','{sample.ts_dna_multiqc}','{sample.ts_pcgr_report}','{sample.ts_pcgr_maf}','{sample.ts_pcgr_snvs_indels}','{sample.ts_pcgr_cna_segments}','{sample.ts_tp_ini}','{sample.ts_run_proc_fastq_1_rna}','{sample.ts_run_proc_fastq_2_rna}','{sample.ts_beluga_fastq_1_rna}','{sample.ts_beluga_fastq_2_rna}','{sample.ts_rna_vcf}','{sample.ts_final_rna_bam_variants}','{sample.ts_rna_multiqc}','{sample.ts_annofuse}','{sample.ts_gridss}','{sample.ts_rna_abundance}','{sample.ts_big_wig_tracks_f}','{sample.ts_big_wig_tracks_r}','{sample.ts_rna_abundance_ini}','{sample.ts_rna_variants_ini}')")

def update_fileloc_details(sample):
    cur = sample.conn.cursor()
    cur.execute(f"""SELECT * FROM File_Locations WHERE Sample='{sample.sample}'""")
    result = cur.fetchone()
    if result:
        cur.execute(f"DELETE FROM File_Locations WHERE Sample='{sample.sample}'")
        cur.execute(f"INSERT INTO File_Locations (Sample,Run_Proc_BAM_DNA_T,Run_Proc_BAM_DNA_N,Beluga_BAM_DNA_T,Beluga_BAM_DNA_N,DNA_VCF_G,DNA_VCF_S,Mutect2_Somatic_vcf,Mutect2_Germline_vcf,strelka2_Germline_vcf,strelka2_Somatic_vcf,vardict_Germline_vcf,vardict_Somatic_vcf,varscan2_Germline_vcf,varscan2_Somatic_vcf,cnvkit_vcf,Final_VCF,Final_DNA_BAM_T,Final_DNA_BAM_N,DNA_MultiQC,pcgr_report,pcgr_maf,pcgr_snvs_indels,pcgr_cna_segments,TP_ini,Run_Proc_fastq_1_RNA,Run_Proc_fastq_2_RNA,Beluga_fastq_1_RNA,Beluga_fastq_2_RNA,RNA_VCF,Final_RNA_BAM_variants,RNA_MultiQC,AnnoFuse,GRIDSS,RNA_Abundance,big_wig_tracks_F,big_wig_tracks_R,RNA_Abundance_ini,RNA_Variants_ini) VALUES ('{sample.sample}','{sample.run_proc_bam_dna_t}','{sample.run_proc_bam_dna_n}','{sample.beluga_bam_dna_t}','{sample.beluga_bam_dna_n}','{sample.dna_vcf_g}','{sample.dna_vcf_s}','{sample.mutect2_somatic_vcf}','{sample.mutect2_germline_vcf}','{sample.strelka2_germline_vcf}','{sample.strelka2_somatic_vcf}','{sample.vardict_germline_vcf}','{sample.vardict_somatic_vcf}','{sample.varscan2_germline_vcf}','{sample.varscan2_somatic_vcf}','{sample.cnvkit_vcf}','{sample.final_vcf}','{sample.final_dna_bam_t}','{sample.final_dna_bam_n}','{sample.dna_multiqc}','{sample.pcgr_report}','{sample.pcgr_maf}','{sample.pcgr_snvs_indels}','{sample.pcgr_cna_segments}','{sample.tp_ini}','{sample.run_proc_fastq_1_rna}','{sample.run_proc_fastq_2_rna}','{sample.beluga_fastq_1_rna}','{sample.beluga_fastq_2_rna}','{sample.rna_vcf}','{sample.final_rna_bam_variants}','{sample.rna_multiqc}','{sample.annofuse}','{sample.gridss}','{sample.rna_abundance}','{sample.big_wig_tracks_f}','{sample.big_wig_tracks_r}','{sample.rna_abundance_ini}','{sample.rna_variants_ini}')")
    else:
        cur.execute(f"INSERT INTO File_Locations (Sample,Run_Proc_BAM_DNA_T,Run_Proc_BAM_DNA_N,Beluga_BAM_DNA_T,Beluga_BAM_DNA_N,DNA_VCF_G,DNA_VCF_S,Mutect2_Somatic_vcf,Mutect2_Germline_vcf,strelka2_Germline_vcf,strelka2_Somatic_vcf,vardict_Germline_vcf,vardict_Somatic_vcf,varscan2_Germline_vcf,varscan2_Somatic_vcf,cnvkit_vcf,Final_VCF,Final_DNA_BAM_T,Final_DNA_BAM_N,DNA_MultiQC,pcgr_report,pcgr_maf,pcgr_snvs_indels,pcgr_cna_segments,TP_ini,Run_Proc_fastq_1_RNA,Run_Proc_fastq_2_RNA,Beluga_fastq_1_RNA,Beluga_fastq_2_RNA,RNA_VCF,Final_RNA_BAM_variants,RNA_MultiQC,AnnoFuse,GRIDSS,RNA_Abundance,big_wig_tracks_F,big_wig_tracks_R,RNA_Abundance_ini,RNA_Variants_ini) VALUES ('{sample.sample}','{sample.run_proc_bam_dna_t}','{sample.run_proc_bam_dna_n}','{sample.beluga_bam_dna_t}','{sample.beluga_bam_dna_n}','{sample.dna_vcf_g}','{sample.dna_vcf_s}','{sample.mutect2_somatic_vcf}','{sample.mutect2_germline_vcf}','{sample.strelka2_germline_vcf}','{sample.strelka2_somatic_vcf}','{sample.vardict_germline_vcf}','{sample.vardict_somatic_vcf}','{sample.varscan2_germline_vcf}','{sample.varscan2_somatic_vcf}','{sample.cnvkit_vcf}','{sample.final_vcf}','{sample.final_dna_bam_t}','{sample.final_dna_bam_n}','{sample.dna_multiqc}','{sample.pcgr_report}','{sample.pcgr_maf}','{sample.pcgr_snvs_indels}','{sample.pcgr_cna_segments}','{sample.tp_ini}','{sample.run_proc_fastq_1_rna}','{sample.run_proc_fastq_2_rna}','{sample.beluga_fastq_1_rna}','{sample.beluga_fastq_2_rna}','{sample.rna_vcf}','{sample.final_rna_bam_variants}','{sample.rna_multiqc}','{sample.annofuse}','{sample.gridss}','{sample.rna_abundance}','{sample.big_wig_tracks_f}','{sample.big_wig_tracks_r}','{sample.rna_abundance_ini}','{sample.rna_variants_ini}')")


# def update_status_db(
#     conn,
#     sample,
#     DNA_N_Transfered,
#     DNA_T_Transfered,
#     Alignment,
#     Variants,
#     Reports,
#     Tumour_Pair_Complete,
#     RNA_Transfrered,
#     rna_abundance,
#     RNA_Pseudoalign_expression,
#     RNA_Alignment_Variant,
#     RNA_Reports,
#     RNA_Complete
#     ):
#     cur = conn.cursor()
#     cur.execute(f"SELECT * FROM STATUS WHERE Sample='{sample}'")
#     result = cur.fetchone()
#     if result:
#         cur.execute(f"DELETE FROM STATUS WHERE Sample='{sample}'")
#         cur.execute(f"INSERT INTO STATUS (Sample,DNA_N_Transfered,DNA_T_Transfered,Alignment,Variants,Reports,Tumour_Pair_Complete,RNA_Transfrered,RNA_Abundance,RNA_Pseudoalign_expression,RNA_Alignment_Variant,RNA_Reports,RNA_Complete) VALUES ('{sample}','{DNA_N_Transfered}','{DNA_T_Transfered}','{Alignment}','{Variants}','{Reports}','{Tumour_Pair_Complete}','{RNA_Transfrered}','{rna_abundance}','{RNA_Pseudoalign_expression}','{RNA_Alignment_Variant}','{RNA_Reports}','{RNA_Complete}')")
#     else:
#         cur.execute(f"INSERT INTO STATUS (Sample,DNA_N_Transfered,DNA_T_Transfered,Alignment,Variants,Reports,Tumour_Pair_Complete,RNA_Transfrered,RNA_Abundance,RNA_Pseudoalign_expression,RNA_Alignment_Variant,RNA_Reports,RNA_Complete) VALUES ('{sample}','{DNA_N_Transfered}','{DNA_T_Transfered}','{Alignment}','{Variants}','{Reports}','{Tumour_Pair_Complete}','{RNA_Transfrered}','{rna_abundance}','{RNA_Pseudoalign_expression}','{RNA_Alignment_Variant}','{RNA_Reports}','{RNA_Complete}')")

def update_status_db(
    conn,
    patient,
    dna_n_transferred,
    dna_t_transferred,
    dna_alignment,
    dna_variant_call,
    dna_report,
    dna_pipeline_execution,
    dna_delivered,
    rna_transferred,
    rna_pipeline_light_execution,
    rna_light_delivered,
    rna_pseudoalignment,
    rna_variant_call,
    rna_report,
    rna_pipeline_execution,
    rna_delivered
    ):
    cur = conn.cursor()
    cur.execute(f"SELECT * FROM status WHERE patient='{patient}'")
    result = cur.fetchone()
    if result:
        cur.execute("""
            UPDATE status
            SET dna_n_transferred = :dna_n_transferred,
            dna_t_transferred = :dna_t_transferred,
            dna_alignment = :dna_alignment,
            dna_variant_call = :dna_variant_call,
            dna_report = :dna_report,
            dna_pipeline_execution = :dna_pipeline_execution,
            dna_delivered = :dna_delivered,
            rna_transferred = :rna_transferred,
            rna_pipeline_light_execution = :rna_pipeline_light_execution,
            rna_light_delivered = :rna_light_delivered,
            rna_pseudoalignment = :rna_pseudoalignment,
            rna_variant_call = :rna_variant_call,
            rna_report = :rna_report,
            rna_pipeline_execution = :rna_pipeline_execution,
            rna_delivered = :rna_delivered
            WHERE
                patient = :patient
            ;
            """,
            {
            "dna_n_transferred": dna_n_transferred,
            "dna_t_transferred": dna_t_transferred,
            "dna_alignment": dna_alignment,
            "dna_variant_call": dna_variant_call,
            "dna_report": dna_report,
            "dna_pipeline_execution": dna_pipeline_execution,
            "dna_delivered": dna_delivered,
            "rna_transferred": rna_transferred,
            "rna_pipeline_light_execution":rna_pipeline_light_execution,
            "rna_light_delivered": rna_light_delivered,
            "rna_pseudoalignment": rna_pseudoalignment,
            "rna_variant_call": rna_variant_call,
            "rna_report": rna_report,
            "rna_pipeline_execution": rna_pipeline_execution,
            "rna_delivered": rna_delivered,
            "patient": patient
            }
            )
    else:
        cur.execute("""
            INSERT INTO status (
                patient,
                dna_n_transferred,
                dna_t_transferred,
                dna_alignment,
                dna_variant_call,
                dna_report,
                dna_pipeline_execution,
                dna_delivered,
                rna_transferred,
                rna_pipeline_light_execution,
                rna_light_delivered,
                rna_pseudoalignment,
                rna_variant_call,
                rna_report,
                rna_pipeline_execution,
                rna_delivered
                )
            VALUES (
                :patient,
                :dna_n_transferred,
                :dna_t_transferred,
                :dna_alignment,
                :dna_variant_call,
                :dna_report,
                :dna_pipeline_execution,
                :dna_delivered,
                :rna_transferred,
                :rna_pipeline_light_execution,
                :rna_light_delivered,
                :rna_pseudoalignment,
                :rna_variant_call,
                :rna_report,
                :rna_pipeline_execution,
                :rna_delivered,
                )
            ;
            """,
            {
            "dna_n_transferred": dna_n_transferred,
            "dna_t_transferred": dna_t_transferred,
            "dna_alignment": dna_alignment,
            "dna_variant_call": dna_variant_call,
            "dna_report": dna_report,
            "dna_pipeline_execution": dna_pipeline_execution,
            "dna_delivered": dna_delivered,
            "rna_transferred": rna_transferred,
            "rna_pipeline_light_execution":rna_pipeline_light_execution,
            "rna_light_delivered":rna_light_delivered,
            "rna_pseudoalignment": rna_pseudoalignment,
            "rna_variant_call": rna_variant_call,
            "rna_report": rna_report,
            "rna_pipeline_execution": rna_pipeline_execution,
            "rna_delivered": rna_delivered,
            "patient": patient
            }
            )
