#!/usr/bin/env python3

import datetime
import os
import sqlite3
import glob
import getpass
import argparse
import re
import typing
from sqlite3 import Error
from sqlalchemy import Column, ForeignKey, Integer, Boolean, String, JSON, DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relation, sessionmaker


class PasswordPromptAction(argparse.Action):
    def __init__(self,
             option_strings,
             dest=None,
             nargs=0,
             default=None,
             required=False,
             type=None,
             metavar=None,
             help=None):
        super(PasswordPromptAction, self).__init__(
             option_strings=option_strings,
             dest=dest,
             nargs=nargs,
             default=default,
             required=required,
             metavar=metavar,
             type=type,
             help=help)

    def __call__(self, parser, args, values, option_string=None):
        password = getpass.getpass()
        setattr(args, self.dest, password)


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

#Replace the other functions with this one
def extract_value(conn, table, sample, column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM {table} WHERE sample='{sample}'")
    result = cur.fetchone()
    # if not result:
    #     result = "NA"
    cur.row_factory = None
    return result


def extract_sample_metrics(conn, sample, column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM key_metric WHERE sample='{sample}'")
    result = cur.fetchone()
    # if not result:
    #     result = "NA"
    cur.row_factory = None
    return result

def extract_patient_field(conn, patient, column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM patient WHERE patient='{patient}'")
    result = cur.fetchone()
    # if not result:
    #     result = "NA"
    cur.row_factory = None
    return result

def get_run_by_sample(conn, sample):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT run FROM sample WHERE sample='{sample}'")
    result = cur.fetchone()
    # if not result:
    #     result = "NA"
    cur.row_factory = None
    return result

def get_run_by_patient(conn, patient):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT run FROM sample WHERE patient='{patient}'")
    result = cur.fetchone()
    # if not result:
    #     result = "NA"
    cur.row_factory = None
    return result

def extract_fileloc_field(conn, patient, column):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"SELECT {column} FROM file_location WHERE patient='{patient}'")
    result = cur.fetchone()
    # if not result:
    #     result = "NA"
    cur.row_factory = None
    return result

def update_value(conn, table, sample, column, value):
    cur = conn.cursor()
    cur.row_factory = lambda cursor, row: row[0]
    cur.execute(f"""UPDATE {table}
        SET {column} = '{value}'
        WHERE sample={sample} AND {column} IS NULL;""")

def update_key_metric_table(
    conn,
    sample,
    dna_bases_over_q30_percent,
    dna_aligned_reads_count,
    dna_dedup_coverage,
    median_insert_size,
    dna_contamination,
    rna_aligned_reads_count,
    rna_exonic_rate,
    rna_ribosomal_contamination_count,
    dna_concordance,
    dna_tumour_purity,
    flag,
    fail
    ):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM key_metric WHERE sample='{sample}'""")
    result = cur.fetchone()
    if result:
        cur.execute("""
            UPDATE key_metric
            SET dna_bases_over_q30_percent = :dna_bases_over_q30_percent,
            dna_aligned_reads_count = :dna_aligned_reads_count,
            dna_dedup_coverage = :dna_dedup_coverage,
            median_insert_size = :median_insert_size,
            dna_contamination = :dna_contamination,
            rna_exonic_rate = :rna_exonic_rate,
            rna_aligned_reads_count = :rna_aligned_reads_count,
            rna_ribosomal_contamination_count = :rna_ribosomal_contamination_count,
            dna_concordance = :dna_concordance,
            dna_tumour_purity = :dna_tumour_purity,
            flag = :flag,
            fail = :fail
            WHERE
                sample = :sample;
            """,
            {
            "dna_bases_over_q30_percent": dna_bases_over_q30_percent,
            "dna_aligned_reads_count": dna_aligned_reads_count,
            "dna_dedup_coverage": dna_dedup_coverage,
            "median_insert_size": median_insert_size,
            "dna_contamination": dna_contamination,
            "rna_exonic_rate": rna_exonic_rate,
            "rna_aligned_reads_count": rna_aligned_reads_count,
            "rna_ribosomal_contamination_count": rna_ribosomal_contamination_count,
            "dna_concordance": dna_concordance,
            "dna_tumour_purity": dna_tumour_purity,
            "flag": flag,
            "fail": fail,
            "sample": sample
            }
            )
        # cur.execute(f"DELETE FROM key_metric WHERE sample='{sample}'")
        # cur.execute(f"INSERT INTO key_metric (sample,dna_bases_over_q30_percent,dna_aligned_reads_count,dna_raw_coverage,dna_dedup_coverage,median_insert_size,dna_duplication_rate,dna_contamination,rna_raw_reads_count,rna_exonic_rate,rna_aligned_reads_count,rna_ribosomal_contamination_count,dna_concordance,dna_tumour_purity,flag,fail,raw_median_insert_size,raw_mean_insert_size) VALUES ('{sample}','{dna_bases_over_q30_percent}','{dna_aligned_reads_count}','{dna_raw_coverage}','{dna_dedup_coverage}','{median_insert_size}','{dna_duplication_rate}','{dna_contamination}','{rna_raw_reads_count}','{rna_exonic_rate}','{rna_aligned_reads_count}','{rna_ribosomal_contamination_count}','{dna_concordance}','{dna_tumour_purity}','{flag}','{fail}','{raw_median_insert_size}','{raw_mean_insert_size}')")
    else:
        cur.execute("""
            INSERT INTO key_metric (sample,dna_bases_over_q30_percent,dna_aligned_reads_count,dna_raw_coverage,dna_dedup_coverage,median_insert_size,raw_duplication_rate,dna_contamination,raw_reads_count,rna_exonic_rate,rna_aligned_reads_count,rna_ribosomal_contamination_count,dna_concordance,dna_tumour_purity,flag,fail,raw_median_insert_size,raw_mean_insert_size)
            VALUES (:sample,:dna_bases_over_q30_percent,:dna_aligned_reads_count,:dna_raw_coverage,:dna_dedup_coverage,:median_insert_size,:raw_duplication_rate,:dna_contamination,:raw_reads_count,:rna_exonic_rate,:rna_aligned_reads_count,:rna_ribosomal_contamination_count,:dna_concordance,:dna_tumour_purity,:flag,:fail,:raw_median_insert_size,:raw_mean_insert_size)
            """,
            {
            "dna_bases_over_q30_percent": dna_bases_over_q30_percent,
            "dna_aligned_reads_count": dna_aligned_reads_count,
            "raw_mean_coverage": None,
            "dna_dedup_coverage": dna_dedup_coverage,
            "median_insert_size": median_insert_size,
            "raw_duplication_rate": None,
            "dna_contamination": dna_contamination,
            "raw_reads_count": None,
            "rna_exonic_rate": rna_exonic_rate,
            "rna_aligned_reads_count": rna_aligned_reads_count,
            "rna_ribosomal_contamination_count": rna_ribosomal_contamination_count,
            "dna_concordance": dna_concordance,
            "dna_tumour_purity": dna_tumour_purity,
            "flag": flag,
            "fail": fail,
            "raw_median_insert_size": None,
            "raw_mean_insert_size": None,
            "sample": sample
            })

def update_key_metric_table_run_metrics(
    conn,
    sample,
    raw_reads_count,
    raw_mean_coverage,
    raw_median_insert_size,
    raw_mean_insert_size,
    raw_duplication_rate,
    flag,
    fail
    ):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM key_metric WHERE sample='{sample}'""")
    result = cur.fetchone()
    if result:
        cur.execute("""
            UPDATE key_metric
            SET raw_mean_coverage = :raw_mean_coverage,
            raw_duplication_rate = :raw_duplication_rate,
            raw_reads_count = :raw_reads_count,
            flag = :flag,
            fail = :fail,
            raw_median_insert_size = :raw_median_insert_size,
            raw_mean_insert_size = :raw_mean_insert_size
            WHERE
                sample = :sample;
            """,
            {
            "raw_mean_coverage": raw_mean_coverage,
            "raw_duplication_rate": raw_duplication_rate,
            "raw_reads_count": raw_reads_count,
            "flag": flag,
            "fail": fail,
            "raw_median_insert_size": raw_median_insert_size,
            "raw_mean_insert_size": raw_mean_insert_size,
            "sample": sample
            }
            )
        # cur.execute(f"DELETE FROM key_metric WHERE sample='{sample}'")
        # cur.execute(f"INSERT INTO key_metric (sample,dna_bases_over_q30_percent,dna_aligned_reads_count,dna_raw_coverage,dna_dedup_coverage,median_insert_size,dna_duplication_rate,dna_contamination,rna_raw_reads_count,rna_exonic_rate,rna_aligned_reads_count,rna_ribosomal_contamination_count,dna_concordance,dna_tumour_purity,flag,fail,raw_median_insert_size,raw_mean_insert_size) VALUES ('{sample}','{dna_bases_over_q30_percent}','{dna_aligned_reads_count}','{dna_raw_coverage}','{dna_dedup_coverage}','{median_insert_size}','{dna_duplication_rate}','{dna_contamination}','{rna_raw_reads_count}','{rna_exonic_rate}','{rna_aligned_reads_count}','{rna_ribosomal_contamination_count}','{dna_concordance}','{dna_tumour_purity}','{flag}','{fail}','{raw_median_insert_size}','{raw_mean_insert_size}')")
    else:
        cur.execute("""
            INSERT INTO key_metric (sample,dna_bases_over_q30_percent,dna_aligned_reads_count,raw_mean_coverage,dna_dedup_coverage,median_insert_size,raw_duplication_rate,dna_contamination,raw_reads_count,rna_exonic_rate,rna_aligned_reads_count,rna_ribosomal_contamination_count,dna_concordance,dna_tumour_purity,flag,fail,raw_median_insert_size,raw_mean_insert_size)
            VALUES (:sample,:dna_bases_over_q30_percent,:dna_aligned_reads_count,:raw_mean_coverage,:dna_dedup_coverage,:median_insert_size,:raw_duplication_rate,:dna_contamination,:raw_reads_count,:rna_exonic_rate,:rna_aligned_reads_count,:rna_ribosomal_contamination_count,:dna_concordance,:dna_tumour_purity,:flag,:fail,:raw_median_insert_size,:raw_mean_insert_size)
            """,
            {
            "dna_bases_over_q30_percent": None,
            "dna_aligned_reads_count": None,
            "raw_mean_coverage": raw_mean_coverage,
            "dna_dedup_coverage": None,
            "median_insert_size": None,
            "raw_duplication_rate": raw_duplication_rate,
            "dna_contamination": None,
            "raw_reads_count": raw_reads_count,
            "rna_exonic_rate": None,
            "rna_aligned_reads_count": None,
            "rna_ribosomal_contamination_count": None,
            "dna_concordance": None,
            "dna_tumour_purity": None,
            "flag": flag,
            "fail": fail,
            "raw_median_insert_size": raw_median_insert_size,
            "raw_mean_insert_size": raw_mean_insert_size,
            "sample": sample
            })

def update_patient_table(
    conn,
    patient,
    patient_corrected,
    institution,
    cohort,
    dna_n,
    dna_n_corrected,
    dna_t,
    dna_t_corrected,
    rna,
    rna_corrected
    ):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM patient WHERE patient='{patient}'""")
    result = cur.fetchone()
    if result:
        if result[4].startswith('MoH') and dna_n == 'NULL':
            dna_n = result[4]
            dna_n_corrected = result[4]
        if result[6].startswith('MoH') and dna_t == 'NULL':
            dna_t = result[6]
            dna_t_corrected = result[6]
        if result[8].startswith('MoH') and rna == 'NULL':
            rna = result[8]
            rna_corrected = result[8]
        cur.execute("""
            UPDATE patient
            SET patient_corrected = :patient_corrected,
                institution = :institution
                cohort = :cohort
                dna_n = :dna_n
                dna_n_corrected = :dna_n_corrected
                dna_t = :dna_t
                dna_t_corrected = :dna_t_corrected
                rna = :rna
                rna_corrected = :rna_corrected
            WHERE
                patient = :patient;
            """,
            {
            "patient_corrected": {patient_corrected},
            "institution": {institution},
            "cohort": {cohort},
            "dna_n": {dna_n},
            "dna_n_corrected": {dna_n_corrected},
            "dna_t": {dna_t},
            "dna_t_corrected": {dna_t_corrected},
            "rna": {rna},
            "rna_corrected": {rna_corrected},
            "patient": {patient}
            }
            )
        # cur.execute(f"DELETE FROM patient WHERE patient='{patient}'")
        # cur.execute(f"INSERT INTO patient (patient,patient_corrected,institution,cohort,dna_n,dna_n_corrected,dna_t,dna_t_corrected,rna,rna_corrected) VALUES ('{patient}','{patient_corrected}','{institution}','{cohort}','{dna_n}','{dna_n_corrected}','{dna_t}','{dna_t_corrected}','{rna}','{rna_corrected}')")
    else:
        cur.execute("""
            INSERT INTO patient (sample,patient_corrected,institution,cohort,dna_n,dna_n_corrected,dna_t,dna_t_corrected,rna,rna_corrected)
            VALUES (:patient,:patient_corrected,:institution,:cohort,:dna_n,:dna_n_corrected,:dna_t,:dna_t_corrected,:rna,:rna_corrected)
            """,
            {
            "patient_corrected": {patient_corrected},
            "institution": {institution},
            "cohort": {cohort},
            "dna_n": {dna_n},
            "dna_n_corrected": {dna_n_corrected},
            "dna_t": {dna_t},
            "dna_t_corrected": {dna_t_corrected},
            "rna": {rna},
            "rna_corrected": {rna_corrected},
            "patient": {patient}
            })

def update_sample_table(
    conn,
    sample,
    patient,
    run
    ):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM sample WHERE sample='{sample}'""")
    result = cur.fetchone()
    if result:
        cur.execute(f"""
            UPDATE sample
            SET patient = '{patient}',
                run = '{run}'
            WHERE
                sample = '{sample}';
            """)
    else:
        cur.execute(f"""
            INSERT INTO sample (sample,patient,run)
            VALUES ('{sample}','{patient}','{run}')
            """)

def update_readset_table(
    conn,
    sample,
    readset,
    library_type,
    run_type,
    run,
    lane,
    adapter1,
    adapter2,
    quality_offset,
    bed,
    fastq1,
    fastq2,
    bam
    ):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM readset WHERE sample='{sample}'""")
    result = cur.fetchone()
    if result:
        cur.execute("""
            UPDATE readset
            SET readset = :readset,
                library_type = :library_type,
                run_type = :run_type,
                run = :run,
                lane = :lane,
                adapter1 = :adapter1,
                adapter2 = :adapter2,
                quality_offset = :quality_offset,
                bed = :bed,
                fastq1 = :fastq1,
                fastq2 = :fastq2,
                bam = :bam
            WHERE
                sample = :sample;
            """,
            {
            "readset": readset,
            "library_type": library_type,
            "run_type": run_type,
            "run": run,
            "lane": lane,
            "adapter1": adapter1,
            "adapter2": adapter2,
            "quality_offset": quality_offset,
            "bed": bed,
            "fastq1": fastq1,
            "fastq2": fastq2,
            "bam": bam,
            "sample": sample
            }
            )
    else:
        cur.execute("""
            INSERT INTO readset (sample,readset,library_type,run_type,run,lane,adapter1,adapter2,quality_offset,bed,fastq1,fastq2,bam)
            VALUES (:sample,:readset,:library_type,:run_type,:run,:lane,:adapter1,:adapter2,:quality_offset,:bed,:fastq1,:fastq2,:bam);
            """,
            {
            "readset": readset,
            "library_type": library_type,
            "run_type": run_type,
            "run": run,
            "lane": lane,
            "adapter1": adapter1,
            "adapter2": adapter2,
            "quality_offset": quality_offset,
            "bed": bed,
            "fastq1": fastq1,
            "fastq2": fastq2,
            "bam": bam,
            "sample": sample
            }
            )

#Grabs all fields from sample table and outputs an array of the data
def extract_patient_details(conn, patient_corrected):
    # Create cursor object
    cur = conn.cursor()
    #Test to see if it exists and if so extract the current data
    cur.execute(f"""SELECT * FROM patient WHERE patient='{patient_corrected}'""")
    result = cur.fetchone()
    # if not result:
        #return an aray that gets turned into an object for all the data from a sample.
        # result = 'NA'
    return result

def extract_patient_names(conn):
    # Create cursor object
    cur = conn.cursor()
    cur.execute("""SELECT patient FROM patient""")
    cur.row_factory = lambda cursor, row: row[0]
    result = cur.fetchall()
    cur.row_factory = None
    if not result:
        result = []
    return result

def extract_sample_names(conn):
    # Create cursor object
    cur = conn.cursor()
    cur.execute("""SELECT sample FROM key_metric""")
    cur.row_factory = lambda cursor, row: row[0]
    result = cur.fetchall()
    cur.row_factory = None
    if not result:
        result = []
    return result

def extract_fileloc_details(conn, patient):
    cur = conn.cursor()
    cur.execute(f"""SELECT * FROM file_location WHERE patient='{patient}'""")
    result = cur.fetchone()
    # if not result:
    #     result = ['NA'] * 50
    return result


def extract_timestamp_details(conn, patient):
    cur = conn.cursor()
    cur.execute(f"""SELECT * FROM timestamp WHERE patient='{patient}'""")
    result = cur.fetchone()
    # if not result:
    #     result = ['NA'] * 50
    return result

def extract_readset_details(conn, sample):
    cur = conn.cursor()
    cur.execute(f"""SELECT * FROM readset WHERE sample='{sample}'""")
    result = cur.fetchone()
    # if not result:
    #     result = ['NA'] * 50
    return result

def update_timestamp_details(patient):
    cur = patient.conn.cursor()
    cur.execute(f"""SELECT * FROM timestamp WHERE patient='{patient.patient}'""")
    result = cur.fetchone()
    if result:
        cur.execute("""
            UPDATE timestamp
            SET dna_abacus_raw_bam_t = :timestamp_dna_abacus_raw_bam_t,
            dna_abacus_raw_bam_n = :timestamp_dna_abacus_raw_bam_n,
            dna_raw_bam_t = :timestamp_dna_raw_bam_t,
            dna_raw_bam_n = :timestamp_dna_raw_bam_n,
            dna_germline_vcf = :timestamp_dna_germline_vcf,
            dna_somatic_vcf = :timestamp_dna_somatic_vcf,
            dna_mutect2_somatic_vcf = :timestamp_dna_mutect2_somatic_vcf,
            dna_mutect2_germline_vcf = :timestamp_dna_mutect2_germline_vcf,
            dna_strelka2_germline_vcf = :timestamp_dna_strelka2_germline_vcf,
            dna_strelka2_somatic_vcf = :timestamp_dna_strelka2_somatic_vcf,
            dna_vardict_germline_vcf = :timestamp_dna_vardict_germline_vcf,
            dna_vardict_somatic_vcf = :timestamp_dna_vardict_somatic_vcf,
            dna_varscan2_germline_vcf = :timestamp_dna_varscan2_germline_vcf,
            dna_varscan2_somatic_vcf = :timestamp_dna_varscan2_somatic_vcf,
            dna_final_vcf = :timestamp_dna_final_vcf,
            dna_final_bam_t = :timestamp_dna_final_bam_t,
            dna_final_bam_n = :timestamp_dna_final_bam_n,
            dna_multiqc_report = :timestamp_dna_multiqc_report,
            dna_pcgr_report = :timestamp_dna_pcgr_report,
            dna_tumour_pair_ini = :timestamp_dna_tumour_pair_ini,
            rna_abacus_raw_fastq1 = :timestamp_rna_abacus_raw_fastq1,
            rna_abacus_raw_fastq2 = :timestamp_rna_abacus_raw_fastq2,
            rna_raw_fastq1 = :timestamp_rna_raw_fastq1,
            rna_raw_fastq2 = :timestamp_rna_raw_fastq2,
            rna_vcf = :timestamp_rna_vcf,
            rna_final_bam_expression = :timestamp_rna_final_bam_expression,
            rna_final_bam_variants = :timestamp_rna_final_bam_variants,
            rna_multiqc_report = :timestamp_rna_multiqc_report,
            rna_annofuse_tsv = :timestamp_rna_annofuse_tsv,
            rna_gridss_report = :timestamp_rna_gridss_report,
            rna_abundance_stringtie = :timestamp_rna_abundance_stringtie,
            rna_forward_bigwig = :timestamp_rna_forward_bigwig,
            rna_reverse_bigwig = :timestamp_rna_reverse_bigwig,
            rna_abundance_ini = :timestamp_rna_abundance_ini,
            rna_variants_ini = :timestamp_rna_variants_ini
            WHERE
                patient = :patient;
            """,
            {
            "dna_abacus_raw_bam_t": patient.timestamp_dna_abacus_raw_bam_t,
            "dna_abacus_raw_bam_n": patient.timestamp_dna_abacus_raw_bam_n,
            "dna_raw_bam_t": patient.timestamp_dna_raw_bam_t,
            "dna_raw_bam_n": patient.timestamp_dna_raw_bam_n,
            "dna_germline_vcf": patient.timestamp_dna_germline_vcf,
            "dna_somatic_vcf": patient.timestamp_dna_somatic_vcf,
            "dna_mutect2_somatic_vcf": patient.timestamp_dna_mutect2_somatic_vcf,
            "dna_mutect2_germline_vcf": patient.timestamp_dna_mutect2_germline_vcf,
            "dna_strelka2_germline_vcf": patient.timestamp_dna_strelka2_germline_vcf,
            "dna_strelka2_somatic_vcf": patient.timestamp_dna_strelka2_somatic_vcf,
            "dna_vardict_germline_vcf": patient.timestamp_dna_vardict_germline_vcf,
            "dna_vardict_somatic_vcf": patient.timestamp_dna_vardict_somatic_vcf,
            "dna_varscan2_germline_vcf": patient.timestamp_dna_varscan2_germline_vcf,
            "dna_varscan2_somatic_vcf": patient.timestamp_dna_varscan2_somatic_vcf,
            "dna_final_vcf": patient.timestamp_dna_final_vcf,
            "dna_final_bam_t": patient.timestamp_dna_final_bam_t,
            "dna_final_bam_n": patient.timestamp_dna_final_bam_n,
            "dna_multiqc_report": patient.timestamp_dna_multiqc_report,
            "dna_pcgr_report": patient.timestamp_dna_pcgr_report,
            "dna_tumour_pair_ini": patient.timestamp_dna_tumour_pair_ini,
            "rna_abacus_raw_fastq1": patient.timestamp_rna_abacus_raw_fastq1,
            "rna_abacus_raw_fastq2": patient.timestamp_rna_abacus_raw_fastq2,
            "rna_raw_fastq1": patient.timestamp_rna_raw_fastq1,
            "rna_raw_fastq2": patient.timestamp_rna_raw_fastq2,
            "rna_vcf": patient.timestamp_rna_vcf,
            "rna_final_bam_expression": patient.timestamp_rna_final_bam_expression,
            "rna_final_bam_variants": patient.timestamp_rna_final_bam_variants,
            "rna_multiqc_report": patient.timestamp_rna_multiqc_report,
            "rna_annofuse_tsv": patient.timestamp_rna_annofuse_tsv,
            "rna_gridss_report": patient.timestamp_rna_gridss_report,
            "rna_abundance_stringtie": patient.timestamp_rna_abundance_stringtie,
            "rna_forward_bigwig": patient.timestamp_rna_forward_bigwig,
            "rna_reverse_bigwig": patient.timestamp_rna_reverse_bigwig,
            "rna_abundance_ini": patient.timestamp_rna_abundance_ini,
            "rna_variants_ini": patient.timestamp_rna_variants_ini,
            "patient": patient.patient
            }
            )
        # cur.execute(f"DELETE FROM timestamp WHERE patient='{patient.patient}'")
        # cur.execute(f"INSERT INTO timestamp (patient,dna_abacus_raw_bam_t,dna_abacus_raw_bam_n,dna_raw_bam_t,dna_raw_bam_n,dna_germline_vcf,dna_somatic_vcf,dna_mutect2_somatic_vcf,dna_mutect2_germline_vcf,dna_strelka2_germline_vcf,dna_strelka2_somatic_vcf,dna_vardict_germline_vcf,dna_vardict_somatic_vcf,dna_varscan2_germline_vcf,dna_varscan2_somatic_vcf,dna_final_vcf,dna_final_bam_t,dna_final_bam_n,dna_multiqc_report,dna_pcgr_report,dna_tumour_pair_ini,rna_abacus_raw_fastq1,rna_abacus_raw_fastq2,rna_raw_fastq1,rna_raw_fastq2,rna_vcf,rna_final_bam_expression,rna_final_bam_variants,rna_multiqc_report,rna_annofuse_tsv,rna_gridss_report,rna_abundance_stringtie,rna_forward_bigwig,rna_reverse_bigwig,rna_abundance_ini,rna_variants_ini) VALUES ('{patient.patient}','{patient.timestamp_dna_abacus_raw_bam_t}','{patient.timestamp_dna_abacus_raw_bam_n}','{patient.timestamp_dna_raw_bam_t}','{patient.timestamp_dna_raw_bam_n}','{patient.timestamp_dna_germline_vcf}','{patient.timestamp_dna_somatic_vcf}','{patient.timestamp_dna_mutect2_somatic_vcf}','{patient.timestamp_dna_mutect2_germline_vcf}','{patient.timestamp_dna_strelka2_germline_vcf}','{patient.timestamp_dna_strelka2_somatic_vcf}','{patient.timestamp_dna_vardict_germline_vcf}','{patient.timestamp_dna_vardict_somatic_vcf}','{patient.timestamp_dna_varscan2_germline_vcf}','{patient.timestamp_dna_varscan2_somatic_vcf}','{patient.timestamp_dna_final_vcf}','{patient.timestamp_dna_final_bam_t}','{patient.timestamp_dna_final_bam_n}','{patient.timestamp_dna_multiqc_report}','{patient.timestamp_dna_pcgr_report}','{patient.timestamp_dna_tumour_pair_ini}','{patient.timestamp_rna_abacus_raw_fastq1}','{patient.timestamp_rna_abacus_raw_fastq2}','{patient.timestamp_rna_raw_fastq1}','{patient.timestamp_rna_raw_fastq2}','{patient.timestamp_rna_vcf}','{patient.timestamp_rna_final_bam_expression}','{patient.timestamp_rna_final_bam_variants}','{patient.timestamp_rna_multiqc_report}','{patient.timestamp_rna_annofuse_tsv}','{patient.timestamp_rna_gridss_report}','{patient.timestamp_rna_abundance_stringtie}','{patient.timestamp_rna_forward_bigwig}','{patient.timestamp_rna_reverse_bigwig}','{patient.timestamp_rna_abundance_ini}','{patient.timestamp_rna_variants_ini}')")
    else:
        cur.execute("""
            INSERT INTO timestamp (patient,dna_abacus_raw_bam_t,dna_abacus_raw_bam_n,dna_raw_bam_t,dna_raw_bam_n,dna_germline_vcf,dna_somatic_vcf,dna_mutect2_somatic_vcf,dna_mutect2_germline_vcf,dna_strelka2_germline_vcf,dna_strelka2_somatic_vcf,dna_vardict_germline_vcf,dna_vardict_somatic_vcf,dna_varscan2_germline_vcf,dna_varscan2_somatic_vcf,dna_final_vcf,dna_final_bam_t,dna_final_bam_n,dna_multiqc_report,dna_pcgr_report,dna_tumour_pair_ini,rna_abacus_raw_fastq1,rna_abacus_raw_fastq2,rna_raw_fastq1,rna_raw_fastq2,rna_vcf,rna_final_bam_expression,rna_final_bam_variants,rna_multiqc_report,rna_annofuse_tsv,rna_gridss_report,rna_abundance_stringtie,rna_forward_bigwig,rna_reverse_bigwig,rna_abundance_ini,rna_variants_ini)
            VALUES (:patient,:timestamp_dna_abacus_raw_bam_t,:timestamp_dna_abacus_raw_bam_n,:timestamp_dna_raw_bam_t,:timestamp_dna_raw_bam_n,:timestamp_dna_germline_vcf,:timestamp_dna_somatic_vcf,:timestamp_dna_mutect2_somatic_vcf,:timestamp_dna_mutect2_germline_vcf,:timestamp_dna_strelka2_germline_vcf,:timestamp_dna_strelka2_somatic_vcf,:timestamp_dna_vardict_germline_vcf,:timestamp_dna_vardict_somatic_vcf,:timestamp_dna_varscan2_germline_vcf,:timestamp_dna_varscan2_somatic_vcf,:timestamp_dna_final_vcf,:timestamp_dna_final_bam_t,:timestamp_dna_final_bam_n,:timestamp_dna_multiqc_report,:timestamp_dna_pcgr_report,:timestamp_dna_tumour_pair_ini,:timestamp_rna_abacus_raw_fastq1,:timestamp_rna_abacus_raw_fastq2,:timestamp_rna_raw_fastq1,:timestamp_rna_raw_fastq2,:timestamp_rna_vcf,:timestamp_rna_final_bam_expression,:timestamp_rna_final_bam_variants,:timestamp_rna_multiqc_report,:timestamp_rna_annofuse_tsv,:timestamp_rna_gridss_report,:timestamp_rna_abundance_stringtie,:timestamp_rna_forward_bigwig,:timestamp_rna_reverse_bigwig,:timestamp_rna_abundance_ini,:timestamp_rna_variants_ini)
            """,
            {
            "dna_abacus_raw_bam_t": patient.timestamp_dna_abacus_raw_bam_t,
            "dna_abacus_raw_bam_n": patient.timestamp_dna_abacus_raw_bam_n,
            "dna_raw_bam_t": patient.timestamp_dna_raw_bam_t,
            "dna_raw_bam_n": patient.timestamp_dna_raw_bam_n,
            "dna_germline_vcf": patient.timestamp_dna_germline_vcf,
            "dna_somatic_vcf": patient.timestamp_dna_somatic_vcf,
            "dna_mutect2_somatic_vcf": patient.timestamp_dna_mutect2_somatic_vcf,
            "dna_mutect2_germline_vcf": patient.timestamp_dna_mutect2_germline_vcf,
            "dna_strelka2_germline_vcf": patient.timestamp_dna_strelka2_germline_vcf,
            "dna_strelka2_somatic_vcf": patient.timestamp_dna_strelka2_somatic_vcf,
            "dna_vardict_germline_vcf": patient.timestamp_dna_vardict_germline_vcf,
            "dna_vardict_somatic_vcf": patient.timestamp_dna_vardict_somatic_vcf,
            "dna_varscan2_germline_vcf": patient.timestamp_dna_varscan2_germline_vcf,
            "dna_varscan2_somatic_vcf": patient.timestamp_dna_varscan2_somatic_vcf,
            "dna_final_vcf": patient.timestamp_dna_final_vcf,
            "dna_final_bam_t": patient.timestamp_dna_final_bam_t,
            "dna_final_bam_n": patient.timestamp_dna_final_bam_n,
            "dna_multiqc_report": patient.timestamp_dna_multiqc_report,
            "dna_pcgr_report": patient.timestamp_dna_pcgr_report,
            "dna_tumour_pair_ini": patient.timestamp_dna_tumour_pair_ini,
            "rna_abacus_raw_fastq1": patient.timestamp_rna_abacus_raw_fastq1,
            "rna_abacus_raw_fastq2": patient.timestamp_rna_abacus_raw_fastq2,
            "rna_raw_fastq1": patient.timestamp_rna_raw_fastq1,
            "rna_raw_fastq2": patient.timestamp_rna_raw_fastq2,
            "rna_vcf": patient.timestamp_rna_vcf,
            "rna_final_bam_expression": patient.timestamp_rna_final_bam_expression,
            "rna_final_bam_variants": patient.timestamp_rna_final_bam_variants,
            "rna_multiqc_report": patient.timestamp_rna_multiqc_report,
            "rna_annofuse_tsv": patient.timestamp_rna_annofuse_tsv,
            "rna_gridss_report": patient.timestamp_rna_gridss_report,
            "rna_abundance_stringtie": patient.timestamp_rna_abundance_stringtie,
            "rna_forward_bigwig": patient.timestamp_rna_forward_bigwig,
            "rna_reverse_bigwig": patient.timestamp_rna_reverse_bigwig,
            "rna_abundance_ini": patient.timestamp_rna_abundance_ini,
            "rna_variants_ini": patient.timestamp_rna_variants_ini,
            "patient": patient.patient
            })

def update_fileloc_details(patient):
    cur = patient.conn.cursor()
    cur.execute(f"""SELECT * FROM file_location WHERE patient='{patient.patient}'""")
    result = cur.fetchone()
    if result:
        cur.execute("""
            UPDATE file_location
            SET dna_abacus_raw_bam_t = :dna_abacus_raw_bam_t,
            dna_abacus_raw_bam_n = :dna_abacus_raw_bam_n,
            dna_raw_bam_t = :dna_raw_bam_t,
            dna_raw_bam_n = :dna_raw_bam_n,
            dna_germline_vcf = :dna_germline_vcf,
            dna_somatic_vcf = :dna_somatic_vcf,
            dna_mutect2_somatic_vcf = :dna_mutect2_somatic_vcf,
            dna_mutect2_germline_vcf = :dna_mutect2_germline_vcf,
            dna_strelka2_germline_vcf = :dna_strelka2_germline_vcf,
            dna_strelka2_somatic_vcf = :dna_strelka2_somatic_vcf,
            dna_vardict_germline_vcf = :dna_vardict_germline_vcf,
            dna_vardict_somatic_vcf = :dna_vardict_somatic_vcf,
            dna_varscan2_germline_vcf = :dna_varscan2_germline_vcf,
            dna_varscan2_somatic_vcf = :dna_varscan2_somatic_vcf,
            dna_final_vcf = :dna_final_vcf,
            dna_final_bam_t = :dna_final_bam_t,
            dna_final_bam_n = :dna_final_bam_n,
            dna_multiqc_report = :dna_multiqc_report,
            dna_pcgr_report = :dna_pcgr_report,
            dna_tumour_pair_ini = :dna_tumour_pair_ini,
            rna_abacus_raw_fastq1 = :rna_abacus_raw_fastq1,
            rna_abacus_raw_fastq2 = :rna_abacus_raw_fastq2,
            rna_raw_fastq1 = :rna_raw_fastq1,
            rna_raw_fastq2 = :rna_raw_fastq2,
            rna_vcf = :rna_vcf,
            rna_final_bam_expression = :rna_final_bam_expression,
            rna_final_bam_variants = :rna_final_bam_variants,
            rna_multiqc_report = :rna_multiqc_report,
            rna_annofuse_tsv = :rna_annofuse_tsv,
            rna_gridss_report = :rna_gridss_report,
            rna_abundance_stringtie = :rna_abundance_stringtie,
            rna_forward_bigwig = :rna_forward_bigwig,
            rna_reverse_bigwig = :rna_reverse_bigwig,
            rna_abundance_ini = :rna_abundance_ini,
            rna_variants_ini = :rna_variants_ini
            WHERE
                patient = :patient;
            """,
            {
            "dna_abacus_raw_bam_t": patient.dna_abacus_raw_bam_t,
            "dna_abacus_raw_bam_n": patient.dna_abacus_raw_bam_n,
            "dna_raw_bam_t": patient.dna_raw_bam_t,
            "dna_raw_bam_n": patient.dna_raw_bam_n,
            "dna_germline_vcf": patient.dna_germline_vcf,
            "dna_somatic_vcf": patient.dna_somatic_vcf,
            "dna_mutect2_somatic_vcf": patient.dna_mutect2_somatic_vcf,
            "dna_mutect2_germline_vcf": patient.dna_mutect2_germline_vcf,
            "dna_strelka2_germline_vcf": patient.dna_strelka2_germline_vcf,
            "dna_strelka2_somatic_vcf": patient.dna_strelka2_somatic_vcf,
            "dna_vardict_germline_vcf": patient.dna_vardict_germline_vcf,
            "dna_vardict_somatic_vcf": patient.dna_vardict_somatic_vcf,
            "dna_varscan2_germline_vcf": patient.dna_varscan2_germline_vcf,
            "dna_varscan2_somatic_vcf": patient.dna_varscan2_somatic_vcf,
            "dna_final_vcf": patient.dna_final_vcf,
            "dna_final_bam_t": patient.dna_final_bam_t,
            "dna_final_bam_n": patient.dna_final_bam_n,
            "dna_multiqc_report": patient.dna_multiqc_report,
            "dna_pcgr_report": patient.dna_pcgr_report,
            "dna_tumour_pair_ini": patient.dna_tumour_pair_ini,
            "rna_abacus_raw_fastq1": patient.rna_abacus_raw_fastq1,
            "rna_abacus_raw_fastq2": patient.rna_abacus_raw_fastq2,
            "rna_raw_fastq1": patient.rna_raw_fastq1,
            "rna_raw_fastq2": patient.rna_raw_fastq2,
            "rna_vcf": patient.rna_vcf,
            "rna_final_bam_expression": patient.rna_final_bam_expression,
            "rna_final_bam_variants": patient.rna_final_bam_variants,
            "rna_multiqc_report": patient.rna_multiqc_report,
            "rna_annofuse_tsv": patient.rna_annofuse_tsv,
            "rna_gridss_report": patient.rna_gridss_report,
            "rna_abundance_stringtie": patient.rna_abundance_stringtie,
            "rna_forward_bigwig": patient.rna_forward_bigwig,
            "rna_reverse_bigwig": patient.rna_reverse_bigwig,
            "rna_abundance_ini": patient.rna_abundance_ini,
            "rna_variants_ini": patient.rna_variants_ini,
            "patient": patient.patient
            }
            )
        # cur.execute(f"DELETE FROM file_location WHERE patient='{patient.patient}'")
        # cur.execute(f"""INSERT INTO file_location (patient,dna_abacus_raw_bam_t,dna_abacus_raw_bam_n,dna_raw_bam_t,dna_raw_bam_n,dna_germline_vcf,dna_somatic_vcf,dna_mutect2_somatic_vcf,dna_mutect2_germline_vcf,dna_strelka2_germline_vcf,dna_strelka2_somatic_vcf,dna_vardict_germline_vcf,dna_vardict_somatic_vcf,dna_varscan2_germline_vcf,dna_varscan2_somatic_vcf,dna_final_vcf,dna_final_bam_t,dna_final_bam_n,dna_multiqc_report,dna_pcgr_report,dna_tumour_pair_ini,rna_abacus_raw_fastq1,rna_abacus_raw_fastq2,rna_raw_fastq1,rna_raw_fastq2,rna_vcf,rna_final_bam_expression,rna_final_bam_variants,rna_multiqc_report,rna_annofuse_tsv,rna_gridss_report,rna_abundance_stringtie,rna_forward_bigwig,rna_reverse_bigwig,rna_abundance_ini,rna_variants_ini) VALUES ('{patient.patient}','{patient.dna_abacus_raw_bam_t}','{patient.dna_abacus_raw_bam_n}','{patient.dna_raw_bam_t}','{patient.dna_raw_bam_n}','{patient.dna_germline_vcf}','{patient.dna_somatic_vcf}','{patient.dna_mutect2_somatic_vcf}','{patient.dna_mutect2_germline_vcf}','{patient.dna_strelka2_germline_vcf}','{patient.dna_strelka2_somatic_vcf}','{patient.dna_vardict_germline_vcf}','{patient.dna_vardict_somatic_vcf}','{patient.dna_varscan2_germline_vcf}','{patient.dna_varscan2_somatic_vcf}','{patient.dna_final_vcf}','{patient.dna_final_bam_t}','{patient.dna_final_bam_n}','{patient.dna_multiqc_report}','{patient.dna_pcgr_report}','{patient.dna_tumour_pair_ini}','{patient.rna_abacus_raw_fastq1}','{patient.rna_abacus_raw_fastq2}','{patient.rna_raw_fastq1}','{patient.rna_raw_fastq2}','{patient.rna_vcf}','{patient.rna_final_bam_expression}','{patient.rna_final_bam_variants}','{patient.rna_multiqc_report}','{patient.rna_annofuse_tsv}','{patient.rna_gridss_report}','{patient.rna_abundance_stringtie}','{patient.rna_forward_bigwig}','{patient.rna_reverse_bigwig}','{patient.rna_abundance_ini}','{patient.rna_variants_ini}')""")
    else:
        cur.execute("""
            INSERT INTO file_location (patient,dna_abacus_raw_bam_t,dna_abacus_raw_bam_n,dna_raw_bam_t,dna_raw_bam_n,dna_germline_vcf,dna_somatic_vcf,dna_mutect2_somatic_vcf,dna_mutect2_germline_vcf,dna_strelka2_germline_vcf,dna_strelka2_somatic_vcf,dna_vardict_germline_vcf,dna_vardict_somatic_vcf,dna_varscan2_germline_vcf,dna_varscan2_somatic_vcf,dna_final_vcf,dna_final_bam_t,dna_final_bam_n,dna_multiqc_report,dna_pcgr_report,dna_tumour_pair_ini,rna_abacus_raw_fastq1,rna_abacus_raw_fastq2,rna_raw_fastq1,rna_raw_fastq2,rna_vcf,rna_final_bam_expression,rna_final_bam_variants,rna_multiqc_report,rna_annofuse_tsv,rna_gridss_report,rna_abundance_stringtie,rna_forward_bigwig,rna_reverse_bigwig,rna_abundance_ini,rna_variants_ini)
            VALUES (:patient,:dna_abacus_raw_bam_t,:dna_abacus_raw_bam_n,:dna_raw_bam_t,:dna_raw_bam_n,:dna_germline_vcf,:dna_somatic_vcf,:dna_mutect2_somatic_vcf,:dna_mutect2_germline_vcf,:dna_strelka2_germline_vcf,:dna_strelka2_somatic_vcf,:dna_vardict_germline_vcf,:dna_vardict_somatic_vcf,:dna_varscan2_germline_vcf,:dna_varscan2_somatic_vcf,:dna_final_vcf,:dna_final_bam_t,:dna_final_bam_n,:dna_multiqc_report,:dna_pcgr_report,:dna_tumour_pair_ini,:rna_abacus_raw_fastq1,:rna_abacus_raw_fastq2,:rna_raw_fastq1,:rna_raw_fastq2,:rna_vcf,:rna_final_bam_expression,:rna_final_bam_variants,:rna_multiqc_report,:rna_annofuse_tsv,:rna_gridss_report,:rna_abundance_stringtie,:rna_forward_bigwig,:rna_reverse_bigwig,:rna_abundance_ini,:rna_variants_ini)
            """,
            {
            "dna_abacus_raw_bam_t": patient.dna_abacus_raw_bam_t,
            "dna_abacus_raw_bam_n": patient.dna_abacus_raw_bam_n,
            "dna_raw_bam_t": patient.dna_raw_bam_t,
            "dna_raw_bam_n": patient.dna_raw_bam_n,
            "dna_germline_vcf": patient.dna_germline_vcf,
            "dna_somatic_vcf": patient.dna_somatic_vcf,
            "dna_mutect2_somatic_vcf": patient.dna_mutect2_somatic_vcf,
            "dna_mutect2_germline_vcf": patient.dna_mutect2_germline_vcf,
            "dna_strelka2_germline_vcf": patient.dna_strelka2_germline_vcf,
            "dna_strelka2_somatic_vcf": patient.dna_strelka2_somatic_vcf,
            "dna_vardict_germline_vcf": patient.dna_vardict_germline_vcf,
            "dna_vardict_somatic_vcf": patient.dna_vardict_somatic_vcf,
            "dna_varscan2_germline_vcf": patient.dna_varscan2_germline_vcf,
            "dna_varscan2_somatic_vcf": patient.dna_varscan2_somatic_vcf,
            "dna_final_vcf": patient.dna_final_vcf,
            "dna_final_bam_t": patient.dna_final_bam_t,
            "dna_final_bam_n": patient.dna_final_bam_n,
            "dna_multiqc_report": patient.dna_multiqc_report,
            "dna_pcgr_report": patient.dna_pcgr_report,
            "dna_tumour_pair_ini": patient.dna_tumour_pair_ini,
            "rna_abacus_raw_fastq1": patient.rna_abacus_raw_fastq1,
            "rna_abacus_raw_fastq2": patient.rna_abacus_raw_fastq2,
            "rna_raw_fastq1": patient.rna_raw_fastq1,
            "rna_raw_fastq2": patient.rna_raw_fastq2,
            "rna_vcf": patient.rna_vcf,
            "rna_final_bam_expression": patient.rna_final_bam_expression,
            "rna_final_bam_variants": patient.rna_final_bam_variants,
            "rna_multiqc_report": patient.rna_multiqc_report,
            "rna_annofuse_tsv": patient.rna_annofuse_tsv,
            "rna_gridss_report": patient.rna_gridss_report,
            "rna_abundance_stringtie": patient.rna_abundance_stringtie,
            "rna_forward_bigwig": patient.rna_forward_bigwig,
            "rna_reverse_bigwig": patient.rna_reverse_bigwig,
            "rna_abundance_ini": patient.rna_abundance_ini,
            "rna_variants_ini": patient.rna_variants_ini,
            "patient": patient.patient
            }
            )

def update_fileloc_details_run_processing(
    patient,
    dna_abacus_raw_bam_t=None,
    dna_abacus_raw_bam_n=None,
    rna_abacus_raw_fastq1=None,
    rna_abacus_raw_fastq2=None
    ):
    cur = patient.conn.cursor()
    cur.execute(f"""SELECT * FROM file_location WHERE patient='{patient.patient}'""")
    result = cur.fetchone()
    if result:
        cur.execute("""
            UPDATE file_location
            SET dna_abacus_raw_bam_t = :dna_abacus_raw_bam_t,
            dna_abacus_raw_bam_n = :dna_abacus_raw_bam_n,
            rna_abacus_raw_fastq1 = :rna_abacus_raw_fastq1,
            rna_abacus_raw_fastq2 = :rna_abacus_raw_fastq2,
            WHERE
                patient = :patient;
            """,
            {
            "dna_abacus_raw_bam_t": dna_abacus_raw_bam_t,
            "dna_abacus_raw_bam_n": dna_abacus_raw_bam_n,
            "rna_abacus_raw_fastq1": rna_abacus_raw_fastq1,
            "rna_abacus_raw_fastq2": rna_abacus_raw_fastq2,
            "patient": patient
            }
            )
        # cur.execute(f"DELETE FROM file_location WHERE patient='{patient.patient}'")
        # cur.execute(f"""INSERT INTO file_location (patient,dna_abacus_raw_bam_t,dna_abacus_raw_bam_n,dna_raw_bam_t,dna_raw_bam_n,dna_germline_vcf,dna_somatic_vcf,dna_mutect2_somatic_vcf,dna_mutect2_germline_vcf,dna_strelka2_germline_vcf,dna_strelka2_somatic_vcf,dna_vardict_germline_vcf,dna_vardict_somatic_vcf,dna_varscan2_germline_vcf,dna_varscan2_somatic_vcf,dna_final_vcf,dna_final_bam_t,dna_final_bam_n,dna_multiqc_report,dna_pcgr_report,dna_tumour_pair_ini,rna_abacus_raw_fastq1,rna_abacus_raw_fastq2,rna_raw_fastq1,rna_raw_fastq2,rna_vcf,rna_final_bam_expression,rna_final_bam_variants,rna_multiqc_report,rna_annofuse_tsv,rna_gridss_report,rna_abundance_stringtie,rna_forward_bigwig,rna_reverse_bigwig,rna_abundance_ini,rna_variants_ini) VALUES ('{patient.patient}','{patient.dna_abacus_raw_bam_t}','{patient.dna_abacus_raw_bam_n}','{patient.dna_raw_bam_t}','{patient.dna_raw_bam_n}','{patient.dna_germline_vcf}','{patient.dna_somatic_vcf}','{patient.dna_mutect2_somatic_vcf}','{patient.dna_mutect2_germline_vcf}','{patient.dna_strelka2_germline_vcf}','{patient.dna_strelka2_somatic_vcf}','{patient.dna_vardict_germline_vcf}','{patient.dna_vardict_somatic_vcf}','{patient.dna_varscan2_germline_vcf}','{patient.dna_varscan2_somatic_vcf}','{patient.dna_final_vcf}','{patient.dna_final_bam_t}','{patient.dna_final_bam_n}','{patient.dna_multiqc_report}','{patient.dna_pcgr_report}','{patient.dna_tumour_pair_ini}','{patient.rna_abacus_raw_fastq1}','{patient.rna_abacus_raw_fastq2}','{patient.rna_raw_fastq1}','{patient.rna_raw_fastq2}','{patient.rna_vcf}','{patient.rna_final_bam_expression}','{patient.rna_final_bam_variants}','{patient.rna_multiqc_report}','{patient.rna_annofuse_tsv}','{patient.rna_gridss_report}','{patient.rna_abundance_stringtie}','{patient.rna_forward_bigwig}','{patient.rna_reverse_bigwig}','{patient.rna_abundance_ini}','{patient.rna_variants_ini}')""")
    else:
        cur.execute("""
            INSERT INTO file_location (patient,dna_abacus_raw_bam_t,dna_abacus_raw_bam_n,dna_raw_bam_t,dna_raw_bam_n,dna_germline_vcf,dna_somatic_vcf,dna_mutect2_somatic_vcf,dna_mutect2_germline_vcf,dna_strelka2_germline_vcf,dna_strelka2_somatic_vcf,dna_vardict_germline_vcf,dna_vardict_somatic_vcf,dna_varscan2_germline_vcf,dna_varscan2_somatic_vcf,dna_final_vcf,dna_final_bam_t,dna_final_bam_n,dna_multiqc_report,dna_pcgr_report,dna_tumour_pair_ini,rna_abacus_raw_fastq1,rna_abacus_raw_fastq2,rna_raw_fastq1,rna_raw_fastq2,rna_vcf,rna_final_bam_expression,rna_final_bam_variants,rna_multiqc_report,rna_annofuse_tsv,rna_gridss_report,rna_abundance_stringtie,rna_forward_bigwig,rna_reverse_bigwig,rna_abundance_ini,rna_variants_ini)
            VALUES (:patient,:dna_abacus_raw_bam_t,:dna_abacus_raw_bam_n,:dna_raw_bam_t,:dna_raw_bam_n,:dna_germline_vcf,:dna_somatic_vcf,:dna_mutect2_somatic_vcf,:dna_mutect2_germline_vcf,:dna_strelka2_germline_vcf,:dna_strelka2_somatic_vcf,:dna_vardict_germline_vcf,:dna_vardict_somatic_vcf,:dna_varscan2_germline_vcf,:dna_varscan2_somatic_vcf,:dna_final_vcf,:dna_final_bam_t,:dna_final_bam_n,:dna_multiqc_report,:dna_pcgr_report,:dna_tumour_pair_ini,:rna_abacus_raw_fastq1,:rna_abacus_raw_fastq2,:rna_raw_fastq1,:rna_raw_fastq2,:rna_vcf,:rna_final_bam_expression,:rna_final_bam_variants,:rna_multiqc_report,:rna_annofuse_tsv,:rna_gridss_report,:rna_abundance_stringtie,:rna_forward_bigwig,:rna_reverse_bigwig,:rna_abundance_ini,:rna_variants_ini)
            """,
            {
            "dna_abacus_raw_bam_t": dna_abacus_raw_bam_t,
            "dna_abacus_raw_bam_n": dna_abacus_raw_bam_n,
            "dna_raw_bam_t": None,
            "dna_raw_bam_n": None,
            "dna_germline_vcf": None,
            "dna_somatic_vcf": None,
            "dna_mutect2_somatic_vcf": None,
            "dna_mutect2_germline_vcf": None,
            "dna_strelka2_germline_vcf": None,
            "dna_strelka2_somatic_vcf": None,
            "dna_vardict_germline_vcf": None,
            "dna_vardict_somatic_vcf": None,
            "dna_varscan2_germline_vcf": None,
            "dna_varscan2_somatic_vcf": None,
            "dna_final_vcf": None,
            "dna_final_bam_t": None,
            "dna_final_bam_n": None,
            "dna_multiqc_report": None,
            "dna_pcgr_report": None,
            "dna_tumour_pair_ini": None,
            "rna_abacus_raw_fastq1": rna_abacus_raw_fastq1,
            "rna_abacus_raw_fastq2": rna_abacus_raw_fastq2,
            "rna_raw_fastq1": None,
            "rna_raw_fastq2": None,
            "rna_vcf": None,
            "rna_final_bam_expression": None,
            "rna_final_bam_variants": None,
            "rna_multiqc_report": None,
            "rna_annofuse_tsv": None,
            "rna_gridss_report": None,
            "rna_abundance_stringtie": None,
            "rna_forward_bigwig": None,
            "rna_reverse_bigwig": None,
            "rna_abundance_ini": None,
            "rna_variants_ini": None,
            "patient": patient
            }
            )

def update_status_db(
    conn,
    patient,
    dna_n_transferred,
    dna_t_transferred,
    dna_alignment,
    dna_variant_call,
    dna_report,
    dna_pipeline_execution,
    rna_transferred,
    rna_alignment,
    rna_variant_call,
    rna_report,
    rna_pipeline_execution
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
            rna_transferred = :rna_transferred,
            rna_alignment = :rna_alignment,
            rna_variant_call = :rna_variant_call,
            rna_report = :rna_report,
            rna_pipeline_execution = :rna_pipeline_execution
            WHERE
                patient = :patient;
            """,
            {
            "dna_n_transferred": dna_n_transferred,
            "dna_t_transferred": dna_t_transferred,
            "dna_alignment": dna_alignment,
            "dna_variant_call": dna_variant_call,
            "dna_report": dna_report,
            "dna_pipeline_execution": dna_pipeline_execution,
            "rna_transferred": rna_transferred,
            "rna_alignment": rna_alignment,
            "rna_variant_call": rna_variant_call,
            "rna_report": rna_report,
            "rna_pipeline_execution": rna_pipeline_execution,
            "patient": patient
            }
            )
        # cur.execute(f"DELETE FROM status WHERE patient='{patient}'")
        # cur.execute(f"INSERT INTO status (patient,dna_n_transferred,dna_t_transferred,dna_alignment,dna_variant_call,dna_report,dna_pipeline_execution,rna_transferred,rna_alignment,rna_variant_call,rna_report,rna_pipeline_execution) VALUES ('{patient}','{dna_n_transferred}','{dna_t_transferred}','{dna_alignment}','{dna_variant_call}','{dna_report}','{dna_pipeline_execution}','{rna_transferred}','{rna_alignment}','{rna_variant_call}','{rna_report}','{rna_pipeline_execution}')")
    else:
        cur.execute("""
            INSERT INTO status (patient,dna_n_transferred,dna_t_transferred,dna_alignment,dna_variant_call,dna_report,dna_pipeline_execution,rna_transferred,rna_alignment,rna_variant_call,rna_report,rna_pipeline_execution)
            VALUES (:patient,:dna_n_transferred,:dna_t_transferred,:dna_alignment,:dna_variant_call,:dna_report,:dna_pipeline_execution,:rna_transferred,:rna_alignment,:rna_variant_call,:rna_report,:rna_pipeline_execution)
            """,
            {
            "dna_n_transferred": dna_n_transferred,
            "dna_t_transferred": dna_t_transferred,
            "dna_alignment": dna_alignment,
            "dna_variant_call": dna_variant_call,
            "dna_report": dna_report,
            "dna_pipeline_execution": dna_pipeline_execution,
            "rna_transferred": rna_transferred,
            "rna_alignment": rna_alignment,
            "rna_variant_call": rna_variant_call,
            "rna_report": rna_report,
            "rna_pipeline_execution": rna_pipeline_execution,
            "patient": patient
            })

def getime(path):
    date = datetime.datetime.fromtimestamp(os.path.getmtime(path))
    return date.strftime("%Y/%m/%d")

class PatientData:
    def __init__(self, connection, patient, beluga_moh_folder):
        # data = []
        self.conn = connection
        self.beluga_moh_folder = beluga_moh_folder
        data = extract_patient_details(connection, patient)
        self.patient = data[0]
        self.patient_corrected = data[1]
        self.institution = data[2]
        self.cohort = data[3]
        self.dna_n = data[4]
        self.dna_n_corrected = data[5]
        self.dna_t = data[6]
        self.dna_t_corrected = data[7]
        self.rna = data[8]
        self.rna_corrected = data[9]
        # self.run = data[10]
        # Defining daughter classes attributes
        self.dna_tumour_pair_ini = None
        self.timestamp_dna_tumour_pair_ini = None
        self.rna_abundance_ini = None
        self.timestamp_rna_abundance_ini = None
        self.rna_variants_ini = None
        self.timestamp_rna_variants_ini = None
        self.dna_abacus_raw_bam_t = None
        self.dna_abacus_raw_bam_n = None
        self.timestamp_dna_abacus_raw_bam_t = None
        self.timestamp_dna_abacus_raw_bam_n = None
        self.rna_abacus_raw_fastq1 = None
        self.rna_abacus_raw_fastq2 = None
        self.timestamp_rna_abacus_raw_fastq1 = None
        self.timestamp_rna_abacus_raw_fastq2 = None
        self.dna_raw_bam_t = None
        self. dna_raw_bam_n= None
        self.timestamp_dna_raw_bam_t = None
        self.timestamp_dna_raw_bam_n = None
        self.rna_raw_fastq1 = None
        self.rna_raw_fastq2 = None
        self.timestamp_rna_raw_fastq1 = None
        self.timestamp_rna_raw_fastq2 = None
        self.dna_final_bam_t = None
        self.dna_final_bam_n = None
        self.timestamp_dna_final_bam_t= None
        self.timestamp_dna_final_bam_n= None
        self.rna_final_bam_expression = None
        self.rna_final_bam_variants = None
        self.timestamp_rna_final_bam_expression = None
        self.timestamp_rna_final_bam_variants = None
        self.dna_germline_vcf = None
        self.dna_somatic_vcf = None
        self.dna_mutect2_somatic_vcf = None
        self.dna_mutect2_germline_vcf = None
        self.dna_strelka2_germline_vcf = None
        self.dna_strelka2_somatic_vcf = None
        self.dna_vardict_germline_vcf = None
        self.dna_vardict_somatic_vcf = None
        self.dna_varscan2_germline_vcf = None
        self.dna_varscan2_somatic_vcf = None
        self.dna_final_vcf = None
        self.timestamp_dna_germline_vcf = None
        self.timestamp_dna_somatic_vcf = None
        self.timestamp_dna_mutect2_somatic_vcf = None
        self.timestamp_dna_mutect2_germline_vcf = None
        self.timestamp_dna_strelka2_germline_vcf = None
        self.timestamp_dna_strelka2_somatic_vcf = None
        self.timestamp_dna_vardict_germline_vcf = None
        self.timestamp_dna_vardict_somatic_vcf = None
        self.timestamp_dna_varscan2_germline_vcf = None
        self.timestamp_dna_varscan2_somatic_vcf = None
        self.timestamp_dna_final_vcf = None
        self.rna_vcf = None
        self.timestamp_rna_vcf= None
        self.dna_multiqc_report = None
        self.dna_pcgr_report = None
        self.timestamp_dna_multiqc_report = None
        self.timestamp_dna_pcgr_report = None
        self.rna_multiqc_report = None
        self.timestamp_rna_multiqc_report = None
        self.rna_annofuse_tsv = None
        self.rna_gridss_report= None
        self.rna_abundance_stringtie= None
        self.rna_forward_bigwig = None
        self.rna_reverse_bigwig = None
        self.timestamp_rna_annofuse_tsv = None
        self.timestamp_rna_gridss_report= None
        self.timestamp_rna_abundance_stringtie= None
        self.timestamp_rna_forward_bigwig = None
        self.timestamp_rna_reverse_bigwig = None

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

class Progress(PatientData):
    def grab_db_values(self):
        # data = []
        data = extract_fileloc_details(self.conn, self.patient)
        if data:
            self.dna_abacus_raw_bam_t = data[1]
            self.dna_abacus_raw_bam_n = data[2]
            self.dna_raw_bam_t = data[3]
            self.dna_raw_bam_n = data[4]
            self.dna_germline_vcf = data[5]
            self.dna_somatic_vcf = data[6]
            self.dna_mutect2_somatic_vcf = data[7]
            self.dna_mutect2_germline_vcf = data[8]
            self.dna_strelka2_germline_vcf = data[9]
            self.dna_strelka2_somatic_vcf = data[10]
            self.dna_vardict_germline_vcf = data[11]
            self.dna_vardict_somatic_vcf = data[12]
            self.dna_varscan2_germline_vcf = data[13]
            self.dna_varscan2_somatic_vcf = data[14]
            self.dna_final_vcf = data[15]
            self.dna_final_bam_t = data[16]
            self.dna_final_bam_n = data[17]
            self.dna_multiqc_report = data[18]
            self.dna_pcgr_report = data[19]
            self.dna_tumour_pair_ini = data[20]
            self.rna_abacus_raw_fastq1 = data[21]
            self.rna_abacus_raw_fastq2 = data[22]
            self.rna_raw_fastq1 = data[23]
            self.rna_raw_fastq2 = data[24]
            self.rna_vcf = data[25]
            self.rna_final_bam_expression = data[26]
            self.rna_final_bam_variants = data[27]
            self.rna_multiqc_report = data[28]
            self.rna_annofuse_tsv = data[29]
            self.rna_gridss_report = data[30]
            self.rna_abundance_stringtie = data[31]
            self.rna_forward_bigwig = data[32]
            self.rna_reverse_bigwig = data[33]
            self.rna_abundance_ini = data[34]
            self.rna_variants_ini = data[35]
        else:
            self.dna_abacus_raw_bam_t = None
            self.dna_abacus_raw_bam_n = None
            self.dna_raw_bam_t = None
            self.dna_raw_bam_n = None
            self.dna_germline_vcf = None
            self.dna_somatic_vcf = None
            self.dna_mutect2_somatic_vcf = None
            self.dna_mutect2_germline_vcf = None
            self.dna_strelka2_germline_vcf = None
            self.dna_strelka2_somatic_vcf = None
            self.dna_vardict_germline_vcf = None
            self.dna_vardict_somatic_vcf = None
            self.dna_varscan2_germline_vcf = None
            self.dna_varscan2_somatic_vcf = None
            self.dna_final_vcf = None
            self.dna_final_bam_t = None
            self.dna_final_bam_n = None
            self.dna_multiqc_report = None
            self.dna_pcgr_report = None
            self.dna_tumour_pair_ini = None
            self.rna_abacus_raw_fastq1 = None
            self.rna_abacus_raw_fastq2 = None
            self.rna_raw_fastq1 = None
            self.rna_raw_fastq2 = None
            self.rna_vcf = None
            self.rna_final_bam_expression = None
            self.rna_final_bam_variants = None
            self.rna_multiqc_report = None
            self.rna_annofuse_tsv = None
            self.rna_gridss_report = None
            self.rna_abundance_stringtie = None
            self.rna_forward_bigwig = None
            self.rna_reverse_bigwig = None
            self.rna_abundance_ini = None
            self.rna_variants_ini = None
        #get the old timestamps
        data = extract_timestamp_details(self.conn, self.patient)
        if data:
            self.timestamp_dna_abacus_raw_bam_t = data[1]
            self.timestamp_dna_abacus_raw_bam_n = data[2]
            self.timestamp_dna_raw_bam_t = data[3]
            self.timestamp_dna_raw_bam_n = data[4]
            self.timestamp_dna_germline_vcf = data[5]
            self.timestamp_dna_somatic_vcf = data[6]
            self.timestamp_dna_mutect2_somatic_vcf = data[7]
            self.timestamp_dna_mutect2_germline_vcf = data[8]
            self.timestamp_dna_strelka2_germline_vcf = data[9]
            self.timestamp_dna_strelka2_somatic_vcf = data[10]
            self.timestamp_dna_vardict_germline_vcf = data[11]
            self.timestamp_dna_vardict_somatic_vcf = data[12]
            self.timestamp_dna_varscan2_germline_vcf = data[13]
            self.timestamp_dna_varscan2_somatic_vcf = data[14]
            self.timestamp_dna_final_vcf = data[15]
            self.timestamp_dna_final_bam_t = data[16]
            self.timestamp_dna_final_bam_n = data[17]
            self.timestamp_dna_multiqc_report = data[18]
            self.timestamp_dna_pcgr_report = data[19]
            self.timestamp_dna_tumour_pair_ini = data[20]
            self.timestamp_rna_abacus_raw_fastq1 = data[21]
            self.timestamp_rna_abacus_raw_fastq2 = data[22]
            self.timestamp_rna_raw_fastq1 = data[23]
            self.timestamp_rna_raw_fastq2 = data[24]
            self.timestamp_rna_vcf = data[25]
            self.timestamp_rna_final_bam_expression = data[26]
            self.timestamp_rna_final_bam_variants = data[27]
            self.timestamp_rna_multiqc_report = data[28]
            self.timestamp_rna_annofuse_tsv = data[29]
            self.timestamp_rna_gridss_report = data[30]
            self.timestamp_rna_abundance_stringtie = data[31]
            self.timestamp_rna_forward_bigwig = data[32]
            self.timestamp_rna_reverse_bigwig = data[33]
            self.timestamp_rna_abundance_ini = data[34]
            self.timestamp_rna_variants_ini = data[35]
        else:
            self.timestamp_dna_abacus_raw_bam_t = None
            self.timestamp_dna_abacus_raw_bam_n = None
            self.timestamp_dna_raw_bam_t = None
            self.timestamp_dna_raw_bam_n = None
            self.timestamp_dna_germline_vcf = None
            self.timestamp_dna_somatic_vcf = None
            self.timestamp_dna_mutect2_somatic_vcf = None
            self.timestamp_dna_mutect2_germline_vcf = None
            self.timestamp_dna_strelka2_germline_vcf = None
            self.timestamp_dna_strelka2_somatic_vcf = None
            self.timestamp_dna_vardict_germline_vcf = None
            self.timestamp_dna_vardict_somatic_vcf = None
            self.timestamp_dna_varscan2_germline_vcf = None
            self.timestamp_dna_varscan2_somatic_vcf = None
            self.timestamp_dna_final_vcf = None
            self.timestamp_dna_final_bam_t = None
            self.timestamp_dna_final_bam_n = None
            self.timestamp_dna_multiqc_report = None
            self.timestamp_dna_pcgr_report = None
            self.timestamp_dna_tumour_pair_ini = None
            self.timestamp_rna_abacus_raw_fastq1 = None
            self.timestamp_rna_abacus_raw_fastq2 = None
            self.timestamp_rna_raw_fastq1 = None
            self.timestamp_rna_raw_fastq2 = None
            self.timestamp_rna_vcf = None
            self.timestamp_rna_final_bam_expression = None
            self.timestamp_rna_final_bam_variants = None
            self.timestamp_rna_multiqc_report = None
            self.timestamp_rna_annofuse_tsv = None
            self.timestamp_rna_gridss_report = None
            self.timestamp_rna_abundance_stringtie = None
            self.timestamp_rna_forward_bigwig = None
            self.timestamp_rna_reverse_bigwig = None
            self.timestamp_rna_abundance_ini = None
            self.timestamp_rna_variants_ini = None

    def __init__(self, connection, patient, beluga_moh_folder):
        super().__init__(connection, patient, beluga_moh_folder)
        self.grab_db_values()

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

    def gather_dna_ini(self):
#################ONCE DONE JUST STORE THE WHOLE FLIPPEN INI HERE ########################################
        # filename = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/TumorPair.config.trace.ini'
        filename = os.path.join(self.beluga_moh_folder, 'MAIN', 'TumorPair.config.trace.ini')
        self.dna_tumour_pair_ini = filename
        self.timestamp_dna_tumour_pair_ini = getime(filename)
#########################################################################################################
    def gather_rna_ini(self):
#################ONCE DONE JUST STORE THE WHOLE FLIPPEN INI HERE########################################
        # filename = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/RnaSeq.config.trace.ini'
        filename = os.path.join(self.beluga_moh_folder, 'MAIN', 'RnaSeq.config.trace.ini')
        self.rna_abundance_ini = filename
        self.timestamp_rna_abundance_ini = getime(filename)
        # filename = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/RnaSeq.config.trace.ini'
        filename = os.path.join(self.beluga_moh_folder, 'MAIN', 'RnaSeq.config.trace.ini')
        self.rna_variants_ini = filename
        self.timestamp_rna_variants_ini = getime(filename)
#########################################################################################################

    def update_status(self):
        dna_n_transferred = None
        if self.timestamp_dna_abacus_raw_bam_n:
            dna_n_transferred = "Complete"
        dna_t_transferred = None
        if self.timestamp_dna_abacus_raw_bam_t:
            dna_t_transferred = "Complete"
        dna_alignment= None
        if self.dna_final_bam_t:
            dna_alignment = "Complete"
        dna_variant_call = None
        if self.dna_germline_vcf and self.dna_somatic_vcf :
            dna_variant_call = "Complete"
        dna_report = None
        if self.dna_multiqc_report and self.dna_pcgr_report:
            dna_report = "Complete"
        dna_pipeline_execution = None
        if dna_report and dna_alignment and dna_variant_call:
            dna_pipeline_execution = "Complete"
        rna_transferred = None
        if self.rna_raw_fastq1:
            rna_transferred = "Complete"
        rna_alignment = None
        if self.rna_final_bam_expression:
            rna_alignment = "Complete"
        rna_variant_call = None
        if self.rna_final_bam_variants:
            rna_variant_call = "Complete"
        rna_report = None
        if self.rna_annofuse_tsv:
            rna_report = "Complete"
        rna_pipeline_execution = None
        if self.rna_reverse_bigwig and rna_alignment and rna_variant_call and rna_report:
            rna_pipeline_execution = "Complete"
        update_status_db(self.conn,self.patient,dna_n_transferred,dna_t_transferred,dna_alignment,dna_variant_call,dna_report,dna_pipeline_execution,rna_transferred,rna_alignment,rna_variant_call,rna_report,rna_pipeline_execution)

    def gather_bam_loc(self):
        # loc1 = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/raw_reads"
        loc1 = os.path.join(self.beluga_moh_folder, 'raw_reads')
        # loc2 = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/raw_reads"
        loc2 = os.path.join(self.beluga_moh_folder, 'MAIN', 'raw_reads')
        if self.dna_n_corrected:
            path = loc1 + "/*" + self.dna_n + "*/*.bam"
            for filename in glob.glob(path):
                self.dna_raw_bam_n = filename
                self.timestamp_dna_raw_bam_n = getime(filename)
            path = loc1 + "/*" + self.dna_t + "*/*.bam"
            for filename in glob.glob(path):
                self.dna_raw_bam_t = filename
                self.timestamp_dna_raw_bam_t = getime(filename)
            path = loc2 + "/*" + self.dna_n + "*/*.bam"
            for filename in glob.glob(path):
                self.dna_raw_bam_n = filename
                self.timestamp_dna_raw_bam_n = getime(filename)
            path = loc2 + "/*" + self.dna_t + "*/*.bam"
            for filename in glob.glob(path):
                self.dna_raw_bam_t = filename
                self.timestamp_dna_raw_bam_t = getime(filename)
        # if self.rna_corrected == "NULL":
        #     self.rna_raw_fastq1 = "NULL"
        #     self.rna_raw_fastq2 = "NULL"
        #     self.timestamp_rna_raw_fastq1 = "NULL"
        #     self.timestamp_rna_raw_fastq2 = "NULL"
        if self.rna_corrected:
            path = loc1 + "/*" + self.rna + "*/*.fastq.gz"
            for filename in glob.glob(path):
                if "R1.fastq.gz" in filename:
                    self.rna_raw_fastq1 = filename
                    self.timestamp_rna_raw_fastq1 = getime(filename)
                elif "R2.fastq.gz" in filename:
                    self.rna_raw_fastq2 = filename
                    self.timestamp_rna_raw_fastq2 = getime(filename)
            path = loc2 + "/*" + self.rna + "*/*.fastq.gz"
            for filename in glob.glob(path):
                if "R1.fastq.gz" in filename:
                    self.rna_raw_fastq1 = filename
                    self.timestamp_rna_raw_fastq1 = getime(filename)
                elif "R2.fastq.gz" in filename:
                    self.rna_raw_fastq2 = filename
                    self.timestamp_rna_raw_fastq2 = getime(filename)

    def gather_final_bams(self):
        if self.dna_n_corrected:
            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.dna_n, "*N.sorted.dup.recal.bam")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'alignment', self.dna_n, f"{self.dna_n}.sorted.dup.recal.bam")
            for filename in glob.glob(path):
                if self.timestamp_dna_final_bam_n != getime(filename):
                    self.dna_final_bam_n = filename
                    self.timestamp_dna_final_bam_n = getime(filename)
                    self.gather_dna_ini()
            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.dna_t, "*T.sorted.dup.recal.bam")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'alignment', self.dna_t, f"{self.dna_t}.sorted.dup.recal.bam")
            for filename in glob.glob(path):
                if self.timestamp_dna_final_bam_t != getime(filename):
                    self.dna_final_bam_t = filename
                    self.timestamp_dna_final_bam_t = getime(filename)
                    self.gather_dna_ini()
        if self.rna_corrected:
            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, "*/*")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'alignment', self.rna, f"{self.rna}.sorted.mdup.*")
            for filename in glob.glob(path):
                if filename.endswith(".sorted.mdup.cram"):
                    if self.timestamp_rna_final_bam_variants != getime(filename):
                        self.timestamp_rna_final_bam_variants = getime(filename)
                        self.rna_final_bam_expression = filename
                elif filename.endswith(".sorted.mdup.split.realigned.recal.bam"):
                    if self.timestamp_rna_final_bam_expression != getime(filename):
                        self.timestamp_rna_final_bam_expression = getime(filename)
                        self.rna_final_bam_variants = filename

    def gather_vcfs(self):
        if self.dna_n_corrected:
            #####################Not Implemented########################################
            self.timestamp_dna_final_vcf = None
            self.dna_final_vcf = None
            ############################################################################
            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/paireddna_variant_call/ensemble", self.patient, "*/*.vt.annot.vcf.gz")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'pairedVariants', 'ensemble', self.patient, f"{self.patient}.ensemble.*.vt.annot.vcf.gz")
            for filename in glob.glob(path):
                if "germline" in filename:
                    if self.timestamp_dna_germline_vcf != getime(filename):
                        self.timestamp_dna_germline_vcf = getime(filename)
                        self.dna_germline_vcf = filename
                elif "somatic" in filename:
                    if self.timestamp_dna_somatic_vcf != getime(filename):
                        self.timestamp_dna_somatic_vcf = getime(filename)
                        self.dna_somatic_vcf = filename
            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/paireddna_variant_call", self.patient, "*/*.vcf.gz")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'pairedVariants', self.patient, f"{self.patient}.*.vcf.gz")
            for filename in glob.glob(path):
                if ".mutect2.somatic.vt.vcf.gz" in filename:
                    if self.timestamp_dna_mutect2_somatic_vcf != getime(filename):
                        self.timestamp_dna_mutect2_somatic_vcf = getime(filename)
                        self.dna_mutect2_somatic_vcf = filename
                if ".mutect2.vcf.gz" in filename:
                    if self.timestamp_dna_mutect2_germline_vcf != getime(filename):
                        self.timestamp_dna_mutect2_germline_vcf = getime(filename)
                        self.dna_mutect2_germline_vcf = filename
                if ".strelka2.somatic.purple.vcf.gz" in filename:
                    if self.timestamp_dna_strelka2_somatic_vcf != getime(filename):
                        self.timestamp_dna_strelka2_somatic_vcf = getime(filename)
                        self.dna_strelka2_somatic_vcf = filename
                if ".strelka2.germline.vt.vcf.gz" in filename:
                    if self.timestamp_dna_strelka2_germline_vcf != getime(filename):
                        self.timestamp_dna_strelka2_germline_vcf = getime(filename)
                        self.dna_strelka2_germline_vcf = filename
                if ".vardict.germline.vt.vcf.gz" in filename:
                    if self.timestamp_dna_vardict_germline_vcf != getime(filename):
                        self.timestamp_dna_vardict_germline_vcf = getime(filename)
                        self.dna_vardict_germline_vcf= filename
                if ".vardict.somatic.vt.vcf.gz" in filename:
                    if self.timestamp_dna_vardict_somatic_vcf != getime(filename):
                        self.timestamp_dna_vardict_somatic_vcf = getime(filename)
                        self.dna_vardict_somatic_vcf = filename
                if ".varscan2.germline.vt.vcf.gz" in filename:
                    if self.timestamp_dna_varscan2_germline_vcf != getime(filename):
                        self.timestamp_dna_varscan2_germline_vcf = getime(filename)
                        self.dna_varscan2_germline_vcf = filename
                if ".varscan2.somatic.vt.vcf.gz" in filename:
                    if self.timestamp_dna_varscan2_somatic_vcf != getime(filename):
                        self.timestamp_dna_varscan2_somatic_vcf = getime(filename)
                        self.dna_varscan2_somatic_vcf = filename
        if self.rna_corrected:
            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment/" + self.rna + "/*.hc.vcf.gz")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'alignment', self.rna, f"{self.rna}.hc.vcf.gz")
            for filename in glob.glob(path):
                if self.timestamp_rna_vcf != getime(filename):
                    self.rna_vcf = filename
                    self.timestamp_rna_vcf = getime(filename)
                    self.gather_rna_ini()


    def gather_reports(self):
        if self.dna_n_corrected:
            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna", self.patient, "*multiqc.html")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'metrics', 'dna', f"{self.patient}.multiqc.html")
            for filename in glob.glob(path):
                if self.timestamp_dna_multiqc_report != getime(filename):
                    self.timestamp_dna_multiqc_report = getime(filename)
                    self.dna_multiqc_report = filename

            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/paireddna_variant_call/ensemble", self.patient, "*/pcgr*/*.flexdb.html")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'pairedVariants', 'ensemble', self.patient, 'pcgr', f"{self.patient}.pcgr_acmg.grch38.flexdb.html")
            for filename in glob.glob(path):
                if self.timestamp_dna_pcgr_report != getime(filename):
                    self.timestamp_dna_pcgr_report = getime(filename)
                    self.dna_pcgr_report = filename

        if self.rna_corrected:
            #####################Not Implemented########################################
            self.rna_multiqc_report = None
            self.timestamp_rna_multiqc_report = None
            ############################################################################


    def gather_rna_other(self):
        if self.rna_corrected:
            #####################Not Implemented########################################
            self.rna_gridss_report= None
            self.timestamp_rna_gridss_report= None
            ############################################################################
            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/fusion", self.rna, "annoFuse/*.putative_driver_fusions.tsv")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'fusion', self.rna, 'annoFuse', f"{self.rna}.putative_driver_fusions.tsv")
            for filename in glob.glob(path):
                if self.timestamp_rna_annofuse_tsv != getime(filename):
                    self.timestamp_rna_annofuse_tsv = getime(filename)
                    self.rna_annofuse_tsv = filename

            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/stringtie", self.rna, "abundance.tab")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'stringtie', self.rna, 'abundance.tab')
            for filename in glob.glob(path):
                if self.timestamp_rna_abundance_stringtie != getime(filename):
                    self.timestamp_rna_abundance_stringtie = getime(filename)
                    self.rna_abundance_stringtie = filename

            # path = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/tracks/bigWig", self.rna, "*")
            path = os.path.join(self.beluga_moh_folder, 'MAIN', 'tracks', 'bigWig', f'{self.rna}.*.bw')
            for filename in glob.glob(path):
                if "forward" in filename:
                    if self.rna_forward_bigwig != getime(filename):
                        self.timestamp_rna_forward_bigwig = getime(filename)
                        self.rna_forward_bigwig = filename
                if "reverse" in filename:
                    if self.rna_reverse_bigwig != getime(filename):
                        self.timestamp_rna_reverse_bigwig = getime(filename)
                        self.rna_reverse_bigwig = filename

def dna_bases_over_q30_percent_check(dna_bases_over_q30_percent, sample_type, fail, flag):
    try:
        if int(dna_bases_over_q30_percent)<75 and sample_type in ('DN', 'DT'):
            fail.append('dna_bases_over_q30_percent')
        elif int(dna_bases_over_q30_percent)<80 and sample_type in ('DN', 'DT'):
            flag.append('dna_bases_over_q30_percent')
    except TypeError:
        pass

def dna_aligned_reads_count_check(dna_aligned_reads_count, sample_type, fail, flag):
    try:
        if int(dna_aligned_reads_count)<260000000 and sample_type == 'DN':
            fail.append('dna_aligned_reads_count')
        elif int(dna_aligned_reads_count)<660000000 and sample_type == 'DN':
            flag.append('dna_aligned_reads_count')
        elif int(dna_aligned_reads_count)<530000000 and sample_type == 'DT':
            fail.append('dna_aligned_reads_count')
        elif int(dna_aligned_reads_count)<1330000000 and sample_type == 'DT':
            flag.append('dna_aligned_reads_count')
    except TypeError:
        pass

def raw_mean_coverage_check(raw_mean_coverage, sample_type, fail):
    try:
        if float(raw_mean_coverage)<30 and sample_type == 'DN':
            fail.append('raw_mean_coverage')
        elif float(raw_mean_coverage)<80 and sample_type == 'DT':
            fail.append('raw_mean_coverage')
    except TypeError:
        pass

def raw_reads_count_check(raw_reads_count, sample_type, fail, flag):
    try:
        if int(raw_reads_count)<80000000 and sample_type == 'RT':
            fail.append('raw_reads_count')
        elif int(raw_reads_count)<100000000 and sample_type == 'RT':
            flag.append('raw_reads_count')
    except TypeError:
        pass

def raw_duplication_rate_check(raw_duplication_rate, sample_type, fail, flag):
    try:
        if float(raw_duplication_rate)>50 and sample_type in ('DT', 'DN'):
            fail.append('raw_duplication_rate')
        elif float(raw_duplication_rate)>20 and sample_type in ('DT', 'DN'):
            flag.append('raw_duplication_rate')
    except TypeError:
        pass

def median_insert_size_check(median_insert_size, fail, flag):
    try:
        if float(median_insert_size)<300:
            flag.append('median_insert_size')
        elif float(median_insert_size)<150:
            fail.append('median_insert_size')
    except TypeError:
        pass

def dna_contamination_check(dna_contamination, fail):
    try:
        if float(dna_contamination)>5:
            fail.append('dna_contamination')
    except TypeError:
        pass

def dna_concordance_check(dna_concordance, fail):
    try:
        if float(dna_concordance)<99:
            fail.append('dna_concordance')
    except TypeError:
        pass

def dna_tumour_purity_check(dna_tumour_purity, fail):
    try:
        if float(dna_tumour_purity)<30:
            fail.append('dna_tumour_purity')
    except TypeError:
        pass

def rna_exonic_rate_check(rna_exonic_rate, fail, flag):
    try:
        if float(rna_exonic_rate)<0.6:
            fail.append('rna_exonic_rate')
        elif float(rna_exonic_rate)<0.8:
            flag.append('rna_exonic_rate')
    except TypeError:
        pass

def rna_ribosomal_contamination_count_check(rrna_count, rna_aligned_reads_count, fail, flag):
    try:
        rna_ribosomal_contamination_count = int(rrna_count)/int(rna_aligned_reads_count)
        if float(rna_ribosomal_contamination_count)>0.35:
            fail.append('rna_ribosomal_contamination_count')
        elif float(rna_ribosomal_contamination_count)>0.1:
            flag.append('rna_ribosomal_contamination_count')
    except (TypeError, ValueError):
        rna_ribosomal_contamination_count = None
    return rna_ribosomal_contamination_count

def extract_rna_ribosomal(sample):
    try:
        # filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample, 'rnaseqc', sample, sample + '.rRNA_counts.txt')
        filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample, 'rnaseqc', sample, sample + '.rRNA_counts.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            line = lines[0]
            fields = line.split("\t")
            ret = fields[0]
    except FileNotFoundError:
        ret = None
    return ret

def parse_rnaseqc_metrics_tmp(sample):
    try:
        # filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample, 'rnaseqc', sample, sample + '.metrics.tmp.txt')
        filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample, 'rnaseqc', sample, sample + '.metrics.tmp.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            rna_aligned_reads_count = lines[3].split("\t")[0]
            rna_exonic_rate = round(float(lines[7].split("\t")[1])*100, 2)
    except FileNotFoundError:
        rna_aligned_reads_count = None
        rna_exonic_rate = None
    return rna_aligned_reads_count, rna_exonic_rate

def parse_run_metrics(sample, run):
    try:
        # filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics', run + '.align_bwa_mem.csv')
        filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics', run + '.align_bwa_mem.csv')
        with open(filename, 'r', encoding="utf-8") as file:
            for line in file:
                parsed_line = line.split(",")
                if sample == parsed_line[6]:
                    raw_reads_count = parsed_line[12]
                    raw_mean_coverage = parsed_line[41]
                    raw_median_insert_size = parsed_line[37]
                    raw_mean_insert_size = parsed_line[38]
                    raw_duplication_rate = parsed_line[15]
    except FileNotFoundError:
        raw_reads_count = None
        raw_mean_coverage = None
        raw_median_insert_size = None
        raw_mean_insert_size = None
        raw_duplication_rate = None
    return raw_reads_count, raw_mean_coverage, raw_median_insert_size, raw_mean_insert_size, raw_duplication_rate

def extract_purity(sample, patient):
    try:
        # filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, 'purple', sample + '.purple.purity.tsv')
        filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, 'purple', sample + '.purple.purity.tsv')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            line = lines[1]
            fields = line.split("\t")
            ret = float(fields[0])*100
    except FileNotFoundError:
        ret = None
    return ret

#WGS N % Contamination,WGS T % Contamination,
def extract_contamination(patient, sample_type):
    if sample_type in ('DN', 'DT'):
        # filename = "".join(glob.glob(os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics', patient + '*DT.contamination.tsv')))
        filename = "".join(glob.glob(os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics', patient + '*DT.contamination.tsv')))
        try:
            with open(filename, 'r', encoding="utf-8") as file:
                for line in file:
                    if line.startswith('Normal') and sample_type == 'DN':
                        ret = line.split(" ")[-1][:-2]
                    elif line.startswith('Tumor') and sample_type == 'DT':
                        ret = line.split(" ")[-1][:-2]
        except FileNotFoundError:
            ret = None
    else:
        ret = None
    return ret

def extract_concordance(patient, sample_type):
    if sample_type in ('DN', 'DT'):
        # filename = "".join(glob.glob(os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics', patient + '*DT.concordance.tsv')))
        filename = "".join(glob.glob(os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics', patient + '*DT.concordance.tsv')))
        try:
            with open(filename, 'r', encoding="utf-8") as file:
                for line in file:
                    if line.startswith('Concordance'):
                        ret = line.split(" ")[-1][:-2]
        except FileNotFoundError:
            ret = None
    else:
        ret = None
    return ret

def extract_insert_size(sample, patient, sample_type):
    try:
        if sample_type == 'RT':
            # filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample + '.insert_size_metrics')
            filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample + '.insert_size_metrics')
            with open(filename, 'r', encoding="utf-8") as file:
                lines = file.readlines()
                line = lines[7]
                metrics = line.split("\t")
                ret = metrics[0]
        elif sample_type in ('DN', 'DT'):
            # filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', patient + '.multiqc_data', 'multiqc_qualimap_bamqc_genome_results.txt')
            filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', patient + '.multiqc_data', 'multiqc_qualimap_bamqc_genome_results.txt')
            with open(filename, 'r', encoding="utf-8") as file:
                for line in file:
                    parsed_line = line.split("\t")
                    if parsed_line[0] == sample:
                        ret = round(float(parsed_line[7]), 0)
        else:
            ret = None
    except FileNotFoundError:
        ret = None
    return ret

def extract_dedup_coverage(sample):
    try:
        # filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'qualimap', sample, 'genome_results.txt')
        filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'qualimap', sample, 'genome_results.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            line = lines[71]
            #line is      mean coverageData = 304.9902X
            metrics = line.split(" ")
            ret = float(metrics[-1].replace('X', ''))
    except (FileNotFoundError, ValueError):
        ret = None
    return ret

def extract_min_aln_rds(sample, patient):
    try:
        # filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', patient + '.multiqc_data', 'multiqc_general_stats.txt')
        filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', patient + '.multiqc_data', 'multiqc_general_stats.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            for line in file:
                parsed_line = line.split("\t")
                if parsed_line[0] == sample:
                    ret = round(float(parsed_line[3]), 0)
    except FileNotFoundError:
        ret = None
    return ret

def extract_bs_over_q30(sample, sample_type):
    try:
        if sample_type in ('DT', 'DN'):
            # filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'picard_metrics', sample + '.all.metrics.quality_distribution_metrics')
            filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'picard_metrics', sample + '.all.metrics.quality_distribution_metrics')
        elif sample_type == 'RT':
            # filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample + '.quality_distribution_metrics')
            filename = os.path.join('/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample + '.quality_distribution_metrics')
        tester = re.compile('(\d+)\W+(\d+)')
        with open(filename, 'r', encoding="utf-8") as file:
            above_30 = 0
            below_30 = 0
            for line in file:
                if line[:1].isdigit():
                    test = tester.match(line)
                    qual = test.group(1)
                    count = test.group(2)
                    if int(qual) < 30:
                        below_30 += int(count)
                    else :
                        above_30 += int(count)
            percent_abv = round((above_30/(above_30+below_30))*100, 2)
    except FileNotFoundError:
        percent_abv = None
    return percent_abv







# Cf. https://stackoverflow.com/questions/48722835/custom-type-hint-annotation
# T = typing.TypeVar('T')

# class json(typing.Generic[T]):
#     pass

# class timestamp(typing.Generic[T]):
#     pass

Base = declarative_base()

class Project(Base):
    """docstring for Project"""
    __tablename__ = "project"

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    def __init__(self, name=None, deleted=None, extra_metadata=None):
        self.name = name
        self.deleted = deleted
        self.extra_metadata = extra_metadata

    def __repr__(self):
        return f"Project({self.name!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Patient(Base):
    """docstring for Patient"""
    __tablename__ = "patient"

    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey("project.id"))
    name = Column(String, nullable=False, unique=True)
    alias = Column(String, nullable=True)
    cohort = Column(String, nullable=True)
    institution = Column(String, nullable=True)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    project = relation("Project", backref="patient", lazy=False)

    def __init__(self, name=None, alias=None, cohort=None, institution=None, deleted=None, extra_metadata=None):
        self.name = name
        self.alias = alias
        self.cohort = cohort
        self.institution = institution
        self.deleted = deleted
        self.extra_metadata = extra_metadata

    def __repr__(self):
        return f"Patient({self.project!r}, {self.name!r}, {self.alias!r}, {self.cohort!r}, {self.institution!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Run(Base):
    """docstring for Patient"""
    __tablename__ = "run"

    id = Column(Integer, primary_key=True)
    lab_id = Column(String, nullable=True)
    name = Column(String, nullable=False, unique=True)
    date = Column(DateTime, nullable=True)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    def __init__(self, lab_id=None, name=None, date=None, deleted=None, extra_metadata=None):
        self.lab_id = lab_id
        self.name = name
        self.date = date
        self.deleted = deleted
        self.extra_metadata = extra_metadata

    def __repr__(self):
        return f"Run({self.lab_id!r}, {self.name!r}, {self.date!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Sample(Base):
    """docstring for Sample"""
    __tablename__ = "sample"

    id = Column(Integer, primary_key=True)
    patient_id = Column(Integer, ForeignKey("patient.id"))
    name = Column(String, nullable=False, unique=True)
    sequencing_technology = Column(String, nullable=True)
    tumour = Column(Boolean, default=False)
    alias = Column(String, nullable=True)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    patient = relation("Patient", backref="sample", lazy=False)

    def __init__(self, name=None, sequencing_technology=None, tumour=None, alias=None, deleted=None, extra_metadata=None):
        self.name = name
        self.sequencing_technology = sequencing_technology
        self.tumour = tumour
        self.alias = alias
        self.deleted = deleted
        self.extra_metadata = extra_metadata

    def __repr__(self):
        return f"Sample({self.patient!r}, {self.name!r}, {self.sequencing_technology!r}, {self.tumour!r}, {self.alias!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Readset(Base):
    """docstring for Readset"""
    __tablename__ = "readset"

    id = Column(Integer, primary_key=True)
    sample_id = Column(Integer, ForeignKey("sample.id"))
    run_id = Column(Integer, ForeignKey("run.id"))
    name = Column(String, nullable=False, unique=True)
    lane = Column(String, nullable=True)
    adapter1 = Column(String, nullable=True)
    adapter2 = Column(String, nullable=True)
    sequencing_type = Column(String, nullable=True)
    quality_offset = Column(String, nullable=True)
    alias = Column(String, nullable=True)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    sample = relation("Sample", backref="readset", lazy=False)
    run = relation("Run", backref="readset", lazy=False)

    def __init__(self, name=None, lane=None, adapter1=None, adapter2=None, sequencing_type=None, quality_offset=None, alias=None, deleted=None, extra_metadata=None):
        self.name = name
        self.lane = lane
        self.adapter1 = adapter1
        self.adapter2 = adapter2
        self.sequencing_type = sequencing_type
        self.quality_offset = quality_offset
        self.alias = alias
        self.deleted = deleted
        self.extra_metadata = extra_metadata

    def __repr__(self):
        return f"Readset({self.sample!r}, {self.run!r}, {self.name!r}, {self.lane!r}, {self.adapter1!r}, {self.adapter2!r}, {self.sequencing_type!r}, {self.quality_offset!r}, {self.alias!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Step(Base):
    """docstring for Step"""
    __tablename__ = "step"

    id = Column(Integer, primary_key=True)
    sample_id = Column(Integer, ForeignKey("sample.id"))
    readset_id = Column(Integer, ForeignKey("readset.id"))
    name = Column(String, nullable=False)
    status = Column(String, nullable=True)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    sample = relation("Sample", backref="step", lazy=False)
    readset = relation("Readset", backref="step", lazy=False)

    def __init__(self, name=None, status=None, deleted=None, extra_metadata=None):
        self.name = name
        self.status = status
        self.deleted = deleted
        self.extra_metadata = extra_metadata

    def __repr__(self):
        return f"Step({self.sample!r}, {self.readset!r}, {self.name!r}, {self.status!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Job(Base):
    """docstring for Job"""
    __tablename__ = "job"

    id = Column(Integer, primary_key=True)
    step_id = Column(Integer, ForeignKey("step.id"))
    name = Column(String, nullable=False)
    start = Column(DateTime, nullable=True)
    stop = Column(DateTime, nullable=True)
    status = Column(String, nullable=True)
    type = Column(String, nullable=True)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    step = relation("Step", backref="job", lazy=False)

    def __init__(self, name=None, start=None, stop=None, status=None, type=None, deleted=None, extra_metadata=None):
        self.name = name
        self.start = start
        self.stop = stop
        self.status = status
        self.type = type
        self.deleted = deleted
        self.extra_metadata = extra_metadata

    def __repr__(self):
        return f"Job({self.step!r}, {self.name!r}, {self.start!r}, {self.stop!r}, {self.status!r}, {self.type!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Metric(Base):
    """docstring for Metric"""
    __tablename__ = "metric"

    id = Column(Integer, primary_key=True)
    job_id = Column(Integer, ForeignKey("job.id"))
    name = Column(String, nullable=False)
    value = Column(String, nullable=True)
    flag = Column(String, nullable=True)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    job = relation("Job", backref="metric", lazy=False)

    def __init__(self, name=None, value=None, flag=None, deleted=None, extra_metadata=None):
        self.name = name
        self.value = value
        self.flag = flag
        self.deleted = deleted
        self.extra_metadata = extra_metadata

    def __repr__(self):
        return f"Metric({self.job!r}, {self.name!r}, {self.value!r}, {self.flag!r}, {self.deleted!r}, {self.extra_metadata!r})"

class File(Base):
    """docstring for File"""
    __tablename__ = "file"

    id = Column(Integer, primary_key=True)
    job_id = Column(Integer, ForeignKey("job.id"))
    path = Column(String, nullable=True)
    type = Column(String, nullable=True)
    description = Column(String, nullable=True)
    creation = Column(DateTime, nullable=True)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    job = relation("Job", backref="file", lazy=False)

    def __init__(self, path=None, type=None, description=None, creation=None, deleted=None, extra_metadata=None):
        self.path = path
        self.type = type
        self.description = description
        self.creation = creation
        self.deleted = deleted
        self.extra_metadata = extra_metadata

    def __repr__(self):
        return f"File({self.job!r}, {self.path!r}, {self.type!r}, {self.description!r}, {self.creation!r}, {self.deleted!r}, {self.extra_metadata!r})"
