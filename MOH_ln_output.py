#!/usr/bin/env python3

# standard import
import argparse
import csv
import datetime
import glob
import logging
import os

import markdown
import pandas as pd
import progressbar
from pymdownx import emoji

from  DB_OPS import (
    create_connection,
    extract_sample_metrics,
    extract_sample_details,
    extract_fileloc_field,
    extract_sample_names,
    extract_patient_status
    )

WIDGETS = [' [', progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') - ', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']

extensions = [
    'markdown.extensions.tables',
    'pymdownx.magiclink',
    'pymdownx.betterem',
    'pymdownx.tilde',
    'pymdownx.emoji',
    'pymdownx.tasklist',
    'pymdownx.superfences',
    'pymdownx.saneheaders',
    'footnotes'
]

extension_configs = {
    "pymdownx.magiclink": {
        "repo_url_shortener": True,
        "repo_url_shorthand": True,
        "provider": "github",
        "user": "facelessuser",
        "repo": "pymdown-extensions"
    },
    "pymdownx.tilde": {
        "subscript": False
    },
    "pymdownx.emoji": {
        "emoji_index": emoji.emojione,
        "emoji_generator": emoji.to_png_sprite,
    }
}

MOH_MAIN_FOLDER = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN"
RAW_READS_FOLDER = os.path.join(MOH_MAIN_FOLDER, "raw_reads")

progressbar.streams.wrap_stderr()
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(prog='MOH_ln_output.py', description="Hardlinks files matching criterias into /lustre03/project/6007512/C3G/projects/share/MOH for delivery.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--black_list', required=False, help="path/to file for patients to be ignored.")
    group.add_argument('--white_list', required=False, help="path/to file for patients to be delivered no matter thresholds.")
    parser.add_argument('--ignore_alignment', required=False, help="Ignores alignment folder for delivery; aka doesn't delivers it.", action='store_true')
    parser.add_argument('--update_metrics', required=False, help="Forces Key_metrics.csv and Warnings.html files generation.", action='store_true')
    parser.add_argument('--update_methods', required=False, help="Forces Methods.html file generation.", action='store_true')
    parser.add_argument('--update_readme', required=False, help="Forces Readme.html file generation.", action='store_true')
    parser.add_argument('--loglevel', help='Sets logging level', choices=logging._levelToName.values(), default='INFO')
    args = parser.parse_args()

    # connection = create_connection("/scratch/stretenp/moh_test/MOH_analysis.db")
    connection = create_connection("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db")
    patients = extract_sample_names(connection)

    if args.black_list:
        black_list = {}
        with open(args.black_list, "r") as black_list_file:
            reader = csv.reader(black_list_file, delimiter="\t")
            for line in reader:
                if line[0] not in black_list:
                    black_list[line[0]] = [line[1]]
                else:
                    black_list[line[0]].append(line[1])
    elif args.white_list:
        white_list = {}
        with open(args.white_list, "r") as white_list_file:
            reader = csv.reader(white_list_file, delimiter="\t")
            for line in reader:
                if line[0] not in white_list:
                    white_list[line[0]] = [line[1]]
                else:
                    white_list[line[0]].append(line[1])
        patients = white_list.keys()
    with progressbar.ProgressBar(max_value=len(patients), widgets=WIDGETS) as progress:
        for index, patient_name in enumerate(patients, 1):
            patient = PatientData(connection, patient_name)
            dna = False
            rna = False
            if args.white_list:
                if extract_patient_status(patient.conn, patient.sample, "dna_pipeline_execution") == "NA":
                    dna = False
                elif "DNA" in white_list[patient_name]:
                    dna = True
                    if float(extract_sample_metrics(patient.conn, patient.dna_n, "WGS_Dedup_Coverage")) <= 30:
                        logger.warning(f"Metric WGS_Dedup_Coverage for Sample {patient.dna_n} is <= 30")
                    if float(extract_sample_metrics(patient.conn, patient.dna_n, "WGS_Dedup_Coverage")) <= 80:
                        logger.warning(f"Metric WGS_Dedup_Coverage for Sample {patient.dna_t} is <= 80")
                if extract_patient_status(patient.conn, patient.sample, "rna_pipeline_light_execution") == "NA":
                    rna = False
                elif "RNA" in white_list[patient_name]:
                    rna = True
                    if extract_sample_metrics(patient.conn, patient.rna, "Raw_Reads_Count") == "NA":
                        logger.warning(f"Metric Raw_Reads_Count for Sample {patient.rna} is NA")
                    elif float(extract_sample_metrics(patient.conn, patient.rna, "Raw_Reads_Count")) <= 80000000:
                        logger.warning(f"Metric Raw_Reads_Count for Sample {patient.rna} is <= 80.000.000")
            else:
                # Check if samples reach the threashold for delivery
                # Check that DNA_T dedup coverage is over 80
                # Check that DNA_N dedup coverage is over 30
                # Check that processing is complete
                # Check that RNA spots is over 100000000
                if args.black_list:
                    if extract_sample_metrics(patient.conn, patient.dna_n, "WGS_Dedup_Coverage") == "NA" or extract_sample_metrics(patient.conn, patient.dna_t, "WGS_Dedup_Coverage") == "NA" or extract_patient_status(patient.conn, patient.sample, "dna_pipeline_execution") == "NA" or (patient_name in black_list and "DNA" in black_list[patient_name]):
                        dna = False
                    elif float(extract_sample_metrics(patient.conn, patient.dna_n, "WGS_Dedup_Coverage")) > 30 and float(extract_sample_metrics(patient.conn, patient.dna_t, "WGS_Dedup_Coverage")) > 80:
                        dna = True
                    if extract_sample_metrics(patient.conn, patient.rna, "Raw_Reads_Count") == "NA" or extract_patient_status(patient.conn, patient.sample, "rna_pipeline_light_execution") == "NA" or (patient_name in black_list and "RNA" in black_list[patient_name]):
                        rna = False
                    elif float(extract_sample_metrics(patient.conn, patient.rna, "Raw_Reads_Count")) > 80000000:
                        rna = True
                else:
                    if extract_sample_metrics(patient.conn, patient.dna_n, "WGS_Dedup_Coverage") == "NA" or extract_sample_metrics(patient.conn, patient.dna_t, "WGS_Dedup_Coverage") == "NA" or extract_patient_status(patient.conn, patient.sample, "dna_pipeline_execution") == "NA":
                        dna = False
                    elif float(extract_sample_metrics(patient.conn, patient.dna_n, "WGS_Dedup_Coverage")) > 30 and float(extract_sample_metrics(patient.conn, patient.dna_t, "WGS_Dedup_Coverage")) > 80:
                        dna = True
                    if extract_sample_metrics(patient.conn, patient.rna, "Raw_Reads_Count") == "NA" or extract_patient_status(patient.conn, patient.sample, "rna_pipeline_light_execution") == "NA":
                        rna = False
                    elif float(extract_sample_metrics(patient.conn, patient.rna, "Raw_Reads_Count")) > 80000000:
                        rna = True

            # Folders used for Delivery
            # base_folder = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/' # Base Folder
            out_folder = '/lustre03/project/6007512/C3G/projects/share/MOH' # Output Folder
            # out_folder = '/scratch/stretenp/MOH' # Output Folder
            # Contains Warnings.txt Readme.txt Log.txt and all subfolders
            out_folder = os.path.join(out_folder, patient.institution, patient.cohort, patient.sample_true)
            # Contains raw bams and fastqs
            raw_folder = os.path.join(out_folder, "raw_data")
            # Contains all variants the subfolder
            var_folder = os.path.join(out_folder, "variants")
            # Contains all the vcfs from the callers
            cal_folder = os.path.join(var_folder, "caller_vcfs")
            # Contains raw cnv
            raw_cnv_folder = os.path.join(out_folder, "raw_cnv")
            # Contains the analysis bams
            alignment_folder = os.path.join(out_folder, "alignment")
            # Contains the ini files
            param_folder = os.path.join(out_folder, "parameters")
            # Contains the reports
            reports_folder = os.path.join(out_folder, "reports")
            pcgr_folder = os.path.join(reports_folder, "pcgr")
            # Contains expression from Kallisto
            expression_folder = os.path.join(out_folder, "expression")
            # Contains all svariants
            svar_folder = os.path.join(out_folder, "svariants")
            linx_folder = os.path.join(svar_folder, "linx")

            # updated keeps track of if we need to update the metrics and warnings file
            # old_log stores the data in the log file to see if things need updating
            updated = False
            old_log = {}
            # See if the directory is created and if so check for file updates.
            log = os.path.join(out_folder, "log.csv")
            warning_file = os.path.join(out_folder, "Warnings.html")
            readme_file = os.path.join(out_folder, "Readme.html")
            methods_file = os.path.join(out_folder, "Methods.html")
            key_metrics_file = os.path.join(out_folder, "Key_metrics.csv")

            os.makedirs(out_folder, exist_ok=True)

            # Load the log file into a dictionary for checking if updates are necessary.
            try:
                with open(log, "r") as log_in:
                    for line in log_in:
                        line.rstrip()
                        fields = line.split(",")
                        old_log[fields[0]] = fields[1]
            except FileNotFoundError:
                with open(log, "w") as log_file:
                    log_file.write("File,Creation,Delivery,Details\n")

            # Populates dna data
            if dna:
                updated = deliver_dna(
                    args,
                    raw_folder,
                    var_folder,
                    svar_folder,
                    linx_folder,
                    cal_folder,
                    raw_cnv_folder,
                    alignment_folder,
                    reports_folder,
                    pcgr_folder,
                    param_folder,
                    connection,
                    patient,
                    log,
                    updated,
                    old_log,
                    )

            # Populates rna data
            if rna:
                updated = deliver_rna(
                    args,
                    raw_folder,
                    expression_folder,
                    var_folder,
                    alignment_folder,
                    reports_folder,
                    pcgr_folder,
                    param_folder,
                    connection,
                    patient,
                    log,
                    updated,
                    old_log,
                    )

            # Not implemented yet
            # if rna and dna:
            #     os.makedirs(var_folder, exist_ok=True)
            #     final_vcf = extract_fileloc_field(connection, sample.sample, "Final_VCF")
            #     updated = get_link_log(final_vcf, var_folder, f"{sample.sample_true}.vcf.gz", log, updated, old_log)

            if args.update_metrics:
                generate_key_metrics(key_metrics_file, log, updated, old_log, patient, connection)
                generate_warning(patient, connection, warning_file, log, updated, old_log)
            if args.update_methods:
                generate_methods(methods_file, log, updated, old_log)
            if args.update_readme:
                generate_readme(readme_file, patient.sample_true, patient.dna_n, patient.dna_t, patient.rna, log, updated, old_log)

            # If any updates were made, generate methods, warning, key metrics and readme files.
            if updated:
                # Add key metrics table for samples
                generate_key_metrics(key_metrics_file, log, updated, old_log, patient, connection)

                # Add warnings file
                generate_warning(patient, connection, warning_file, log, updated, old_log)

                # Add methods file
                generate_methods(methods_file, log, updated, old_log)

                # Add readme file
                generate_readme(readme_file, patient.sample_true, patient.dna_n, patient.dna_t, patient.rna, log, updated, old_log)

            progress.update(index)


def deliver_dna(
    args,
    raw_folder,
    var_folder,
    svar_folder,
    linx_folder,
    cal_folder,
    raw_cnv_folder,
    alignment_folder,
    reports_folder,
    pcgr_folder,
    param_folder,
    connection,
    patient,
    log,
    updated,
    old_log,
    ):

    os.makedirs(raw_folder, exist_ok=True)
    # DNA_N
    if len(glob.glob(os.path.join(RAW_READS_FOLDER, f"{patient.sample}-*DN", "*.bam"))) == 1 and glob.glob(os.path.join(raw_folder, f"{patient.dna_n}.bam")):
        raw_dna_n_bam = os.path.basename(glob.glob(os.path.join(RAW_READS_FOLDER, f"{patient.sample}-*DN", "*.bam"))[0])
        delivered_dna_n_bam = os.path.basename(glob.glob(os.path.join(raw_folder, f"{patient.dna_n}.bam"))[0])
        if raw_dna_n_bam != delivered_dna_n_bam:
            logger.warning(f"The bam file {raw_dna_n_bam} for the sample {patient.dna_n} is already delivered with a different name: {delivered_dna_n_bam}")
    # check if topup
    for bam_n in glob.glob(os.path.join(RAW_READS_FOLDER, f"{patient.sample}-*DN", "*.bam")):
        updated = get_link_log(bam_n, raw_folder, os.path.basename(bam_n), log, updated, old_log)
        bam_n_index = bam_n + ".bai"
        if os.path.exists(bam_n_index):
            updated = get_link_log(bam_n_index, raw_folder, f"{os.path.basename(bam_n)}.bai", log, updated, old_log)
    # DNA_T
    if len(glob.glob(os.path.join(RAW_READS_FOLDER, f"{patient.sample}-*DT", "*.bam"))) == 1 and glob.glob(os.path.join(raw_folder, f"{patient.dna_t}.bam")):
        raw_dna_t_bam = os.path.basename(glob.glob(os.path.join(RAW_READS_FOLDER, f"{patient.sample}-*DT", "*.bam"))[0])
        delivered_dna_t_bam = os.path.basename(glob.glob(os.path.join(raw_folder, f"{patient.dna_t}.bam"))[0])
        if raw_dna_t_bam != delivered_dna_t_bam:
            logger.warning(f"The bam file {raw_dna_t_bam} for the sample {patient.dna_t} is already delivered with a different name: {delivered_dna_t_bam}")
    # check if topup
    for bam_t in glob.glob(os.path.join(RAW_READS_FOLDER, f"{patient.sample}-*DT", "*.bam")):
        updated = get_link_log(bam_t, raw_folder, os.path.basename(bam_t), log, updated, old_log)
        bam_t_index = bam_t + ".bai"
        if os.path.exists(bam_t_index):
            updated = get_link_log(bam_t_index, raw_folder, f"{os.path.basename(bam_t)}.bai", log, updated, old_log)

    os.makedirs(var_folder, exist_ok=True)
    dna_vcf_g = extract_fileloc_field(connection, patient.sample, "DNA_VCF_G")
    updated = get_link_log(dna_vcf_g, var_folder, f"{patient.sample_true}.ensemble.germline.vt.annot.vcf.gz", log, updated, old_log)
    if dna_vcf_g != "NA":
        dna_vcf_g_index = dna_vcf_g + ".tbi"
        if os.path.exists(dna_vcf_g_index):
            updated = get_link_log(dna_vcf_g_index, var_folder, f"{patient.sample_true}.ensemble.germline.vt.annot.vcf.gz.tbi", log, updated, old_log)
        dna_vcf_g_md5 = dna_vcf_g + ".md5"
        if os.path.exists(dna_vcf_g_md5):
            updated = get_link_log(dna_vcf_g_md5, var_folder, f"{patient.sample_true}.ensemble.germline.vt.annot.vcf.gz.md5", log, updated, old_log)
    dna_vcf_s = extract_fileloc_field(connection, patient.sample, "DNA_VCF_S")
    updated = get_link_log(dna_vcf_s, var_folder, f"{patient.sample_true}.ensemble.somatic.vt.annot.vcf.gz", log, updated, old_log)
    if dna_vcf_s != "NA":
        dna_vcf_s_index = dna_vcf_s + ".tbi"
        if os.path.exists(dna_vcf_s_index):
            updated = get_link_log(dna_vcf_s_index, var_folder, f"{patient.sample_true}.ensemble.somatic.vt.annot.vcf.gz.tbi", log, updated, old_log)
        dna_vcf_s_md5 = dna_vcf_s + ".md5"
        if os.path.exists(dna_vcf_s_md5):
            updated = get_link_log(dna_vcf_s_md5, var_folder, f"{patient.sample_true}.ensemble.somatic.vt.annot.vcf.gz.md5", log, updated, old_log)

    os.makedirs(cal_folder, exist_ok=True)
    mutect2_somatic_vcf = extract_fileloc_field(connection, patient.sample, "Mutect2_Somatic_vcf")
    updated = get_link_log(mutect2_somatic_vcf, cal_folder, f"{patient.sample_true}.mutect2.somatic.vt.vcf.gz", log, updated, old_log)
    strelka2_germline_vcf = extract_fileloc_field(connection, patient.sample, "strelka2_Germline_vcf")
    updated = get_link_log(strelka2_germline_vcf, cal_folder, f"{patient.sample_true}.strelka2.germline.vt.vcf.gz", log, updated, old_log)
    strelka2_somatic_vcf = extract_fileloc_field(connection, patient.sample, "strelka2_Somatic_vcf")
    updated = get_link_log(strelka2_somatic_vcf, cal_folder, f"{patient.sample_true}.strelka2.somatic.vt.vcf.gz", log, updated, old_log)
    vardict_germline_vcf = extract_fileloc_field(connection, patient.sample, "vardict_Germline_vcf")
    updated = get_link_log(vardict_germline_vcf, cal_folder, f"{patient.sample_true}.vardict.germline.vt.vcf.gz", log, updated, old_log)
    vardict_somatic_vcf = extract_fileloc_field(connection, patient.sample, "vardict_Somatic_vcf")
    updated = get_link_log(vardict_somatic_vcf, cal_folder, f"{patient.sample_true}.vardict.somatic.vt.vcf.gz", log, updated, old_log)
    varscan2_germline_vcf = extract_fileloc_field(connection, patient.sample, "varscan2_Germline_vcf")
    updated = get_link_log(varscan2_germline_vcf, cal_folder, f"{patient.sample_true}.varscan2.germline.vt.vcf.gz", log, updated, old_log)
    varscan2_somatic_vcf = extract_fileloc_field(connection, patient.sample, "varscan2_Somatic_vcf")
    updated = get_link_log(varscan2_somatic_vcf, cal_folder, f"{patient.sample_true}.varscan2.somatic.vt.vcf.gz", log, updated, old_log)

    os.makedirs(raw_cnv_folder, exist_ok=True)
    cnvkit_vcf = extract_fileloc_field(connection, patient.sample, "cnvkit_vcf")
    updated = get_link_log(cnvkit_vcf, raw_cnv_folder, f"{patient.sample_true}.cnvkit.vcf.gz", log, updated, old_log)
    if cnvkit_vcf != "NA":
        cnvkit_vcf_index = cnvkit_vcf + ".tbi"
        if os.path.exists(cnvkit_vcf_index):
            updated = get_link_log(cnvkit_vcf_index, raw_cnv_folder, f"{patient.sample_true}.cnvkit.vcf.gz.tbi", log, updated, old_log)

    if not args.ignore_alignment:
        os.makedirs(alignment_folder, exist_ok=True)
        final_dna_bam_n = extract_fileloc_field(connection, patient.sample, "Final_DNA_BAM_N")
        updated = get_link_log(final_dna_bam_n, alignment_folder, f"{patient.dna_n}.bam", log, updated, old_log)
        if final_dna_bam_n != "NA":
            final_dna_bam_n_index = final_dna_bam_n + ".bai"
            if os.path.exists(final_dna_bam_n_index):
                updated = get_link_log(final_dna_bam_n_index, alignment_folder, f"{patient.dna_n}.bam.bai", log, updated, old_log)
            final_dna_bam_n_md5 = final_dna_bam_n + ".md5"
            if os.path.exists(final_dna_bam_n_md5):
                updated = get_link_log(final_dna_bam_n_md5, alignment_folder, f"{patient.dna_n}.bam.md5", log, updated, old_log)
        final_dna_bam_t = extract_fileloc_field(connection, patient.sample, "Final_DNA_BAM_T")
        updated = get_link_log(final_dna_bam_t, alignment_folder, f"{patient.dna_t}.bam", log, updated, old_log)
        if final_dna_bam_t != "NA":
            final_dna_bam_t_index = final_dna_bam_t + ".bai"
            if os.path.exists(final_dna_bam_t_index):
                updated = get_link_log(final_dna_bam_t_index, alignment_folder, f"{patient.dna_t}.bam.bai", log, updated, old_log)
            final_dna_bam_t_md5 = final_dna_bam_t + ".md5"
            if os.path.exists(final_dna_bam_t_md5):
                updated = get_link_log(final_dna_bam_t_md5, alignment_folder, f"{patient.dna_t}.bam.md5", log, updated, old_log)

    os.makedirs(reports_folder, exist_ok=True)
    dna_multiqc = extract_fileloc_field(connection, patient.sample, "DNA_MultiQC")
    updated = get_link_log(dna_multiqc, reports_folder, f"{patient.sample_true}_D.multiqc.html", log, updated, old_log)
    pcgr_report = extract_fileloc_field(connection, patient.sample, "pcgr_report")
    updated = get_link_log(pcgr_report, reports_folder, f"{patient.sample_true}_D.pcgr.html", log, updated, old_log)

    os.makedirs(svar_folder, exist_ok=True)
    gridss = extract_fileloc_field(connection, patient.sample, "GRIDSS")
    updated = get_link_log(gridss, svar_folder, f"{patient.sample_true}.gridss", log, updated, old_log)
    gripss_somatic = extract_fileloc_field(connection, patient.sample, "gripss_somatic")
    updated = get_link_log(gripss_somatic, svar_folder, f"{patient.sample_true}.gripss.filtered.somatic.vcf.gz", log, updated, old_log)
    gripss_germline = extract_fileloc_field(connection, patient.sample, "gripss_germline")
    updated = get_link_log(gripss_germline, svar_folder, f"{patient.sample_true}.gripss.filtered.germline.vcf.gz", log, updated, old_log)
    purple_somatic = extract_fileloc_field(connection, patient.sample, "purple_somatic")
    updated = get_link_log(purple_somatic, svar_folder, f"{patient.sample_true}.driver.catalog.somatic.tsv", log, updated, old_log)
    purple_germline = extract_fileloc_field(connection, patient.sample, "purple_germline")
    updated = get_link_log(purple_germline, svar_folder, f"{patient.sample_true}.driver.catalog.germline.tsv", log, updated, old_log)
    purple_circos = extract_fileloc_field(connection, patient.sample, "purple_circos")
    updated = get_link_log(purple_circos, svar_folder, f"{patient.sample_true}.circos.png", log, updated, old_log)
    # linx all tsvs and pngs
    os.makedirs(linx_folder, exist_ok=True)
    linx_genpipes_folder = os.path.join(MOH_MAIN_FOLDER, "SVariants", f"{patient.sample}", "linx")
    for linx_tsv in glob.glob(os.path.join(linx_genpipes_folder, "*.tsv")):
        updated = get_link_log(linx_tsv, linx_folder, os.path.basename(linx_tsv), log, updated, old_log)
    for linx_png in glob.glob(os.path.join(linx_genpipes_folder, "plot", "*.png")):
        updated = get_link_log(linx_png, linx_folder, os.path.basename(linx_png), log, updated, old_log)

    os.makedirs(pcgr_folder, exist_ok=True)
    pcgr_maf = extract_fileloc_field(connection, patient.sample, "pcgr_maf")
    updated = get_link_log(pcgr_maf, pcgr_folder, f"{patient.sample_true}_D.acmg.grch38.maf", log, updated, old_log)
    pcgr_snvs_indels = extract_fileloc_field(connection, patient.sample, "pcgr_snvs_indels")
    updated = get_link_log(pcgr_snvs_indels, pcgr_folder, f"{patient.sample_true}_D.acmg.grch38.snvs_indels.tiers.tsv", log, updated, old_log)
    pcgr_cna_segments = extract_fileloc_field(connection, patient.sample, "pcgr_cna_segments")
    updated = get_link_log(pcgr_cna_segments, pcgr_folder, f"{patient.sample_true}_D.acmg.grch38.cna_segments.tsv.gz", log, updated, old_log)

    os.makedirs(param_folder, exist_ok=True)
    tp_ini = extract_fileloc_field(connection, patient.sample, "TP_ini")
    updated = get_link_log(tp_ini, param_folder, f"{patient.sample_true}.TumourPair.ini", log, updated, old_log)

    return updated


def deliver_rna(
    args,
    raw_folder,
    expression_folder,
    var_folder,
    alignment_folder,
    reports_folder,
    pcgr_folder,
    param_folder,
    connection,
    patient,
    log,
    updated,
    old_log,
    ):
    os.makedirs(raw_folder, exist_ok=True)
    if len(glob.glob(os.path.join(RAW_READS_FOLDER, f"{patient.sample}-*RT", "*.fastq.gz"))) == 2 and glob.glob(os.path.join(raw_folder, f"{patient.rna}*.fastq.gz")):
        for fastq in glob.glob(os.path.join(RAW_READS_FOLDER, f"{patient.sample}-*RT", "*.fastq.gz")):
            if os.path.basename(fastq) not in [os.path.basename(fastq) for fastq in glob.glob(os.path.join(raw_folder, f"{patient.rna}*.fastq.gz"))]:
                logger.warning(f"The fastq file {os.path.basename(fastq)} for the sample {patient.rna} is already delivered with a different name in {raw_folder}")
    for fastq in glob.glob(os.path.join(RAW_READS_FOLDER, f"{patient.sample}-*RT", "*.fastq.gz")):
        updated = get_link_log(fastq, raw_folder, os.path.basename(fastq), log, updated, old_log)

    os.makedirs(expression_folder, exist_ok=True)
    rna_abundance = extract_fileloc_field(connection, patient.sample, "RNA_Abundance")
    updated = get_link_log(rna_abundance, expression_folder, f"{patient.rna}.abundance_transcripts.tsv", log, updated, old_log)
    if rna_abundance != "NA":
        rna_abundance_genes = rna_abundance.replace("transcripts", "genes")
        if os.path.exists(rna_abundance_genes):
            updated = get_link_log(rna_abundance_genes, expression_folder, f"{patient.rna}.abundance_genes.tsv", log, updated, old_log)

    os.makedirs(var_folder, exist_ok=True)
    rna_vcf = extract_fileloc_field(connection, patient.sample, "RNA_VCF")
    updated = get_link_log(rna_vcf, var_folder, f"{patient.rna}.hc.vt.annot.vcf.gz", log, updated, old_log)
    if rna_vcf != "NA":
        rna_vcf_index = rna_vcf + ".tbi"
        if os.path.exists(rna_vcf_index):
            updated = get_link_log(rna_vcf_index, var_folder, f"{patient.rna}.hc.vt.annot.vcf.gz.tbi", log, updated, old_log)
        rna_vcf_md5 = rna_vcf + ".md5"
        if os.path.exists(rna_vcf_md5):
            updated = get_link_log(rna_vcf_md5, var_folder, f"{patient.rna}.hc.vt.annot.vcf.gz.md5", log, updated, old_log)

    rna_vcf_filt = extract_fileloc_field(connection, patient.sample, "RNA_VCF_filt")
    updated = get_link_log(rna_vcf_filt, var_folder, f"{patient.rna}.hc.vt.annot.filt.vcf.gz", log, updated, old_log)
    if rna_vcf_filt != "NA":
        rna_vcf_filt_index = rna_vcf_filt + ".tbi"
        if os.path.exists(rna_vcf_filt_index):
            updated = get_link_log(rna_vcf_filt_index, var_folder, f"{patient.rna}.hc.vt.annot.filt.vcf.gz.tbi", log, updated, old_log)
        rna_vcf_filt_md5 = rna_vcf_filt + ".md5"
        if os.path.exists(rna_vcf_filt_md5):
            updated = get_link_log(rna_vcf_filt_md5, var_folder, f"{patient.rna}.hc.vt.annot.filt.vcf.gz.md5", log, updated, old_log)

    if not args.ignore_alignment:
        os.makedirs(alignment_folder, exist_ok=True)
        final_rna_bam_variants = extract_fileloc_field(connection, patient.sample, "Final_RNA_BAM")
        updated = get_link_log(final_rna_bam_variants, alignment_folder, f"{patient.rna}.variants.bam", log, updated, old_log)
        if final_rna_bam_variants != "NA":
            final_rna_bam_index = final_rna_bam_variants + ".bai"
            if os.path.exists(final_rna_bam_index):
                updated = get_link_log(final_rna_bam_index, alignment_folder, f"{patient.rna}.variants.bam.bai", log, updated, old_log)
            final_rna_bam_md5 = final_rna_bam_variants + ".md5"
            if os.path.exists(final_rna_bam_md5):
                updated = get_link_log(final_rna_bam_md5, alignment_folder, f"{patient.rna}.variants.bam.md5", log, updated, old_log)

    os.makedirs(reports_folder, exist_ok=True)
    rna_multiqc = extract_fileloc_field(connection, patient.sample, "RNA_MultiQC")
    updated = get_link_log(rna_multiqc, reports_folder, f"{patient.sample_true}_R.multiqc.html", log, updated, old_log)
    annofuse = extract_fileloc_field(connection, patient.sample, "AnnoFuse")
    updated = get_link_log(annofuse, reports_folder, f"{patient.rna}.anno_fuse.tsv", log, updated, old_log)

    os.makedirs(param_folder, exist_ok=True)
    rna_expression_ini = extract_fileloc_field(connection, patient.sample, "RNA_Abundance_ini")
    updated = get_link_log(rna_expression_ini, param_folder, f"{patient.sample_true}.RNA.Light.ini", log, updated, old_log)
    rna_variants_ini = extract_fileloc_field(connection, patient.sample, "RNA_Variants_ini")
    updated = get_link_log(rna_variants_ini, param_folder, f"{patient.sample_true}.RNA.Variants.ini", log, updated, old_log)

    # PCGR RNA
    os.makedirs(reports_folder, exist_ok=True)
    rna_pcgr_report = extract_fileloc_field(connection, patient.sample, "rna_pcgr_report")
    updated = get_link_log(rna_pcgr_report, reports_folder, f"{patient.sample_true}_R.pcgr.html", log, updated, old_log)
    os.makedirs(pcgr_folder, exist_ok=True)
    rna_pcgr_maf = extract_fileloc_field(connection, patient.sample, "rna_pcgr_maf")
    updated = get_link_log(rna_pcgr_maf, pcgr_folder, f"{patient.sample_true}_R.acmg.grch38.maf", log, updated, old_log)
    rna_pcgr_snvs_indels = extract_fileloc_field(connection, patient.sample, "rna_pcgr_snvs_indels")
    updated = get_link_log(rna_pcgr_snvs_indels, pcgr_folder, f"{patient.sample_true}_R.acmg.grch38.snvs_indels.tiers.tsv", log, updated, old_log)

    return updated


# Key function. Basiclly updates the logs and links the files based on the database.
def get_link_log(input_file, output_folder, output_file, log, updated, old_log):
    # data = extract_fileloc_field(connection, name, column)
    if input_file != "NA":
        new_time = getime(input_file)
        # Portion for updating files if they have been modified
        if os.path.exists(os.path.join(output_folder, output_file)):
            try:
                old_time = old_log[os.path.join(os.path.basename(output_folder), output_file)]
            except KeyError:
                old_time = getime(input_file)
                log_new(os.path.join(os.path.basename(output_folder), output_file), log, new_time, "New")
            if old_time != new_time:
                os.remove(os.path.join(output_folder, output_file))
                os.link(input_file, os.path.join(output_folder, output_file))
                # os.symlink(input_file, os.path.join(output_folder, output_file))
                log_new(os.path.join(os.path.basename(output_folder), output_file), log, new_time, "Updated")
                updated = True
        else:
            os.link(input_file, os.path.join(output_folder, output_file))
            # os.symlink(input_file, os.path.join(output_folder, output_file))
            log_new(os.path.join(os.path.basename(output_folder), output_file), log, new_time, "New")
            updated = True
    return updated

def get_local_file_log(file, log, updated, old_log):
    #Portion for updating files if they have been modified
    if os.path.exists(file):
        new_time = getime(file)
        old_time = old_log[os.path.basename(file)]
        if old_time != new_time:
            log_new(os.path.basename(file), log, new_time, "Updated")
            updated = True
    else:
        log_new(os.path.basename(file), log, datetime.datetime.today().strftime("%Y/%m/%d"), "New")
        updated = True
    return updated

def log_new(file, log, file_date, message):
    """Update this!!!"""
    if not file_date:
        file_date = datetime.date.today()
        file_date = file_date.strftime("%Y/%m/%d")
    today = datetime.date.today()
    date = today.strftime("%Y/%m/%d")
    # print(file + "," + file_date + "," + date + "," + f"{message}")
    data = file + "," + file_date + "," + date + "," + f"{message}\n"
    with open(log, "a") as log_file:
        log_file.write(data)

def getime(path):
    """Finding timestamp from file"""
    date = datetime.datetime.fromtimestamp(os.path.getmtime(path))
    return date.strftime("%Y/%m/%d")

class PatientData:
    def __init__(self, connection, patient):
        self.conn = connection
        sample_details = extract_sample_details(connection, patient)
        if sample_details["Sample"] == "NA":
            raise Exception(f'No database entry for {patient}')
        self.sample = sample_details["Sample"]
        self.sample_true = sample_details["Sample_True"]
        self.institution = sample_details["Instituion"]
        self.cohort = sample_details["Cohort"]
        self.dna_n = sample_details["DNA_N"]
        self.dna_n_true = sample_details["DNA_N_True"]
        self.dna_t = sample_details["DNA_T"]
        self.dna_t_true = sample_details["DNA_T_True"]
        self.rna = sample_details["RNA"]
        self.rna_true = sample_details["RNA_True"]
    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

def generate_key_metrics(key_metrics_file, log, updated, old_log, patient, connection):
    get_local_file_log(key_metrics_file, log, updated, old_log)
    metrics = pd.read_sql_query(f'select * from KEY_METRICS where Sample="{patient.dna_n}" or Sample="{patient.dna_t}" or Sample="{patient.rna}"', connection)
    metrics.to_csv(key_metrics_file, index=False)

def generate_warning(patient, connection, warning_file, log, updated, old_log):
    warnings_l = []
    warnings = pd.read_sql_query(f'select Sample,Flags,Fails from KEY_METRICS where Sample="{patient.dna_n}" or Sample="{patient.dna_t}" or Sample="{patient.rna}"', connection)
    for _, row in warnings.iterrows():
        sample_name = row["Sample"]
        flags = " - ".join(row["Flags"].split(";"))
        fails = " - ".join(row["Fails"].split(";"))
        warnings_l.append(f"| {sample_name} | &nbsp; {flags} &nbsp; | &nbsp; {fails} |")
    get_local_file_log(warning_file, log, updated, old_log)
    warnings_l_joined = "\n".join(warnings_l)
    data = f"""Below are three columns, "Flags" indicates values that may be troublesome while "Fails" indicates a point of failure. Data may be useable when marked as "Flags", but "Fails" marked data should be carefully considered. If "NA" is present, this data exceeded all standards. Data will not be delivered when coverage is labelled "Fails" at this time.

|  Sample | &nbsp; Flags &nbsp; | &nbsp; Fails |
| :------ | :------------------ | :----------- |
{warnings_l_joined}
"""
    html = markdown.markdown(data, extensions=extensions, extension_configs=extension_configs)
    with open(warning_file, "w", encoding="utf-8", errors="xmlcharrefreplace") as out_file:
        out_file.write(html)

def file_exist_check(file):
    ret = ""
    if os.path.exists(file):
        ret = ":white_check_mark:"
    else:
        ret = ":clock2:"
    return ret

def generate_readme(readme_file, patient, dna_n, dna_t, rna, log, updated, old_log):
    get_local_file_log(readme_file, log, updated, old_log)
    # Add timestamp
    out_folder = os.path.dirname(readme_file)
    timestamp = datetime.datetime.today().strftime("%Y/%m/%d")
    dna_n_raw = ""
    for dna_n_raw_bam in glob.glob(os.path.join(out_folder, "raw_data", "*DN*.bam")):
        dna_n_raw += f"""\n    * `{os.path.basename(dna_n_raw_bam)}` *Raw DNA reads for the Normal sample* {file_exist_check(os.path.join(out_folder, "raw_data", os.path.basename(dna_n_raw_bam)))}"""
    dna_t_raw = ""
    for dna_t_raw_bam in glob.glob(os.path.join(out_folder, "raw_data", "*DT*.bam")):
        dna_t_raw += f"""\n    * `{os.path.basename(dna_t_raw_bam)}` *Raw DNA reads for the Tumor sample* {file_exist_check(os.path.join(out_folder, "raw_data", os.path.basename(dna_t_raw_bam)))}"""
    rna_raw = ""
    for rna_raw_fq in glob.glob(os.path.join(out_folder, "raw_data", "*RT*.fastq.gz")):
        rna_raw += f"""\n    * `{os.path.basename(rna_raw_fq)}` *Raw RNA reads for the Tumor sample* {file_exist_check(os.path.join(out_folder, "raw_data", os.path.basename(rna_raw_fq)))}"""
    linx = ""
    for linx_file in glob.glob(os.path.join(out_folder, "svariants", "linx", "*")):
        linx += f"""\n        * `{os.path.basename(linx_file)}` {file_exist_check(os.path.join(out_folder, "svariants", "linx", os.path.basename(linx_file)))}"""
    data = f"""This directory contains the delivered data for **{patient}** processed by the Canadian Centre for Computational Genomics.
The data will be updated as it becomes available and as such many files may be missing from RNA or DNA upon initial creation of this directory
Should you have concerns, questions, or suggestions, please contact the analysis team at moh-q@computationalgenomics.ca

Within this directory you will find the results of the analysis for a single patient; if multiple bam files are present in the `raw_data` folder it means TopTup(s) have been necessary (when :white_check_mark: is present the file is available, when :clock2: is present the file is coming soon):

* `Readme.html` *This file* :white_check_mark:
* `log.csv` *Log file containing the dates of transfers and if files have been updated* {file_exist_check(os.path.join(out_folder, "log.csv"))}
* [`Warnings.html`](Warnings.html) *Contains details of any warnings and whether they caused a failure of this analysis* {file_exist_check(os.path.join(out_folder, "Warnings.html"))}
* [`Methods.html`](Methods.html) *Contains details on pipelines and references used for the analysis* {file_exist_check(os.path.join(out_folder, "Methods.html"))}
* `Key_metrics.csv` *File with metrics for the patient in csv format* {file_exist_check(os.path.join(out_folder, "Key_metrics.csv"))}
* `raw_data/` *Contains all of the bam's/fastqs from the sequencer. BAM files here include both mapped and unmapped reads and can be converted to the FASTQ format with tools such as SamToFastq.*{dna_n_raw}{dna_t_raw}{rna_raw}
* `variants/` *Contains the vcfs related to variant calls*
    * `{patient}.ensemble.germline.vt.annot.vcf.gz` *Germline Variants found in any of the callers for DNA Sample* {file_exist_check(os.path.join(out_folder, "variants", f"{patient}.ensemble.germline.vt.annot.vcf.gz"))}
    * `{patient}.ensemble.germline.vt.annot.vcf.gz.tbi` *Index of Germline Variants found in any of the callers for DNA Sample* {file_exist_check(os.path.join(out_folder, "variants", f"{patient}.ensemble.germline.vt.annot.vcf.gz.tbi"))}
    * `{patient}.ensemble.somatic.vt.annot.vcf.gz` *Somatic Variants found in any of the callers for DNA Sample* {file_exist_check(os.path.join(out_folder, "variants", f"{patient}.ensemble.somatic.vt.annot.vcf.gz"))}
    * `{patient}.ensemble.somatic.vt.annot.vcf.gz.tbi` *Index of Somatic Variants found in any of the callers for DNA Sample* {file_exist_check(os.path.join(out_folder, "variants", f"{patient}.ensemble.somatic.vt.annot.vcf.gz.tbi"))}
    * `{patient}.hc.vt.annot.vcf.gz` *Variants found using RNA sample* {file_exist_check(os.path.join(out_folder, "variants", f"{patient}.hc.vt.annot.vcf.gz"))}
    * `{patient}.hc.vt.annot.vcf.gz.tbi` *Index of Variants found using RNA sample* {file_exist_check(os.path.join(out_folder, "variants", f"{patient}.hc.vt.annot.vcf.gz.tbi"))}
    * `{patient}.hc.vt.annot.filt.vcf.gz` *Variants found using RNA sample; annotated with RNAEdits and filtered with coverage >=10x and VAF >=5%* {file_exist_check(os.path.join(out_folder, "variants", f"{patient}.hc.vt.annot.filt.vcf.gz"))}
    * `{patient}.hc.vt.annot.filt.vcf.gz.tbi` *Index of Variants found using RNA sample; annotated with RNAEdits and filtered with coverage >=10x and VAF >=5%* {file_exist_check(os.path.join(out_folder, "variants", f"{patient}.hc.vt.annot.filt.vcf.gz.tbi"))}
    * `caller_vcfs/` *Contains the vcfs produced from individual callers on the DNA samples*
        * `{patient}.mutect2.somatic.vt.vcf.gz` *Somatic results for mutect2* {file_exist_check(os.path.join(out_folder, "variants", "caller_vcfs", f"{patient}.mutect2.somatic.vt.vcf.gz"))}
        * `{patient}.strelka2.germline.vt.vcf.gz` *Germline results for strelka2* {file_exist_check(os.path.join(out_folder, "variants", "caller_vcfs", f"{patient}.strelka2.germline.vt.vcf.gz"))}
        * `{patient}.strelka2.somatic.vt.vcf.gz` *Somatic results for strelka2* {file_exist_check(os.path.join(out_folder, "variants", "caller_vcfs", f"{patient}.strelka2.somatic.vt.vcf.gz"))}
        * `{patient}.vardict.germline.vt.vcf.gz` *Germline results for vardict* {file_exist_check(os.path.join(out_folder, "variants", "caller_vcfs", f"{patient}.vardict.germline.vt.vcf.gz"))}
        * `{patient}.vardict.somatic.vt.vcf.gz` *Somatic results for vardict* {file_exist_check(os.path.join(out_folder, "variants", "caller_vcfs", f"{patient}.vardict.somatic.vt.vcf.gz"))}
        * `{patient}.varscan2.germline.vt.vcf.gz` *Germline results for varscan2* {file_exist_check(os.path.join(out_folder, "variants", "caller_vcfs", f"{patient}.varscan2.germline.vt.vcf.gz"))}
        * `{patient}.varscan2.somatic.vt.vcf.gz` *Somatic results for varscan2* {file_exist_check(os.path.join(out_folder, "variants", "caller_vcfs", f"{patient}.varscan2.somatic.vt.vcf.gz"))}
* `svariants/` *Contains the files related to structural variants*
    * `{patient}.gridss.vcf.gz` *Unfiltered structural variant calls with GRIDSS* {file_exist_check(os.path.join(out_folder, "svariants", f"{patient}.gridss.vcf.gz"))}
    * `{patient}.gripss.filtered.somatic.vcf.gz` *Annotated and filtered somatic structural variant calls from GRIDSS using GRIPSS* {file_exist_check(os.path.join(out_folder, "svariants", f"{patient}.gripss.filtered.somatic.vcf.gz"))}
    * `{patient}.gripss.filtered.germline.vcf.gz` *Annotated and filtered germline structural variant calls from GRIDSS using GRIPSS* {file_exist_check(os.path.join(out_folder, "svariants", f"{patient}.gripss.filtered.germline.vcf.gz"))}
    * `{patient}.driver.catalog.somatic.tsv` *Driver somatic structural variant calls* {file_exist_check(os.path.join(out_folder, "svariants", f"{patient}.driver.catalog.somatic.tsv"))}
    * `{patient}.driver.catalog.germline.tsv` *Driver germline structural variant calls* {file_exist_check(os.path.join(out_folder, "svariants", f"{patient}.driver.catalog.germline.tsv"))}
    * [`{patient}.circos.png`](svariants/{patient}.circos.png) *Circos plot of all variants. Cf. https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md#circos* {file_exist_check(os.path.join(out_folder, "svariants", f"{patient}.circos.png"))}
    * `linx/` *Contains structural variant annotated and visualized by LINX Cf. https://github.com/hartwigmedical/hmftools/tree/master/linx*{linx}
* `raw_cnv/` *Contains the raw copy number calls for each patient DNA*
    * `{patient}.cnvkit.vcf.gz` *Raw cnvkit output* {file_exist_check(os.path.join(out_folder, "raw_cnv", f"{patient}.cnvkit.vcf.gz"))}
    * `{patient}.cnvkit.vcf.gz.tbi` *Index of Raw cnvkit output* {file_exist_check(os.path.join(out_folder, "raw_cnv", f"{patient}.cnvkit.vcf.gz.tbi"))}
* `alignment/` *Contains the alignment data for each sample*
    * `{dna_n}.bam` *Alignment of normal against the reference* {file_exist_check(os.path.join(out_folder, "alignment", f"{dna_n}.bam"))}
    * `{dna_n}.bam.bai` *Index of Alignment of normal against the reference* {file_exist_check(os.path.join(out_folder, "alignment", f"{dna_n}.bam.bai"))}
    * `{dna_t}.bam` *Alignment of tumor against the reference* {file_exist_check(os.path.join(out_folder, "alignment", f"{dna_t}.bam"))}
    * `{dna_t}.bam.bai` *Index of Alignment of tumor against the reference* {file_exist_check(os.path.join(out_folder, "alignment", f"{dna_t}.bam.bai"))}
    * `{rna}.variants.bam` *Alignment of tumor RNA against the reference used in variants analysis* {file_exist_check(os.path.join(out_folder, "alignment", f"{rna}.variants.bam"))}
    * `{rna}.variants.bam.bai` *Index of Alignment of tumor RNA against the reference used in variants analysis* {file_exist_check(os.path.join(out_folder, "alignment", f"{rna}.variants.bam.bai"))}
* `expression/` *Contains the transcripts and genes abundance estimation from Kallisto*
    * `{rna}.abundance_transcripts.tsv` *Table with transcript abundance from Kallisto* {file_exist_check(os.path.join(out_folder, "expression", f"{rna}.abundance_transcripts.tsv"))}
    * `{rna}.abundance_genes.tsv` *Table with gene abundance from Kallisto* {file_exist_check(os.path.join(out_folder, "expression", f"{rna}.abundance_genes.tsv"))}
* `reports/` *Contains the reports for the experiment*
    * [`{patient}_D.multiqc.html`](reports/{patient}_D.multiqc.html) *QC report for the DNA analysis* {file_exist_check(os.path.join(out_folder, "reports", f"{patient}_D.multiqc.html"))}
    * [`{patient}_R.multiqc.html`](reports/{patient}_R.multiqc.html) *QC report for the RNA analysis* {file_exist_check(os.path.join(out_folder, "reports", f"{patient}_R.multiqc.html"))}
    * `{rna}.anno_fuse.tsv` *TSV for fusions detected using RN, to be loaded into Excel sheet for easier display* {file_exist_check(os.path.join(out_folder, "reports", f"{rna}.anno_fuse.tsv"))}
    * [`{patient}_D.pcgr.html`](reports/{patient}_D.pcgr.html) *Personal Cancer Genome Reporter report for the DNA analysis; Cf. https://pcgr.readthedocs.io/en/latest* {file_exist_check(os.path.join(out_folder, "reports", f"{patient}_D.pcgr.html"))}
    * [`{patient}_R.pcgr.html`](reports/{patient}_R.pcgr.html) *Personal Cancer Genome Reporter report for the RNA analysis; Cf. https://pcgr.readthedocs.io/en/latest* {file_exist_check(os.path.join(out_folder, "reports", f"{patient}_R.pcgr.html"))}
    * `pcgr/` *Contains raw tables used to generate PCGR reports*
        * `{patient}_D.acmg.grch38.maf` {file_exist_check(os.path.join(out_folder, "reports", "pcgr", f"{patient}_D.acmg.grch38.maf"))}
        * `{patient}_D.acmg.grch38.snvs_indels.tiers.tsv` {file_exist_check(os.path.join(out_folder, "reports", "pcgr", f"{patient}_D.acmg.grch38.snvs_indels.tiers.tsv"))}
        * `{patient}_D.acmg.grch38.cna_segments.tsv.gz` {file_exist_check(os.path.join(out_folder, "reports", "pcgr", f"{patient}_D.acmg.grch38.cna_segments.tsv.gz"))}
        * `{patient}_R.acmg.grch38.maf` {file_exist_check(os.path.join(out_folder, "reports", "pcgr", f"{patient}_R.acmg.grch38.maf"))}
        * `{patient}_R.acmg.grch38.snvs_indels.tiers.tsv` {file_exist_check(os.path.join(out_folder, "reports", "pcgr", f"{patient}_R.acmg.grch38.snvs_indels.tiers.tsv"))}
* `parameters/` *Contains the records of all the Parameters used in the pipeline analysis*
    * `{patient}.TumourPair.ini` *Parameters used in the tumor pair analysis* {file_exist_check(os.path.join(out_folder, "parameters", f"{patient}.TumourPair.ini"))}
    * `{patient}.RNA.Light.ini` *Parameters used in the RNA expression analysis* {file_exist_check(os.path.join(out_folder,"parameters", f"{patient}.RNA.Light.ini"))}
    * `{patient}.RNA.Variants.ini` *Parameters used in the RNA variant analysis* {file_exist_check(os.path.join(out_folder,"parameters", f"{patient}.RNA.Variants.ini"))}

Generated {timestamp}."""
    html = markdown.markdown(data, extensions=extensions, extension_configs=extension_configs)
    with open(readme_file, "w", encoding="utf-8", errors="xmlcharrefreplace") as out_file:
        out_file.write(html)

def generate_methods(methods_file, log, updated, old_log):
    get_local_file_log(methods_file, log, updated, old_log)
    data = """# Methods
Using an Illumina NovaSeq 6000 instrument, Whole Genome Sequencing (WGS) was performed on the tumor and matched normal samples, with a target depth of coverage of 80X and 30X respectively. Similarly, Whole Transcriptome Sequencing (WTS) was done on tumour samples, with a target of 100 million paired-end reads per sample.

Bioinformatics analyses were performed using the [GenPipes][GenPipes_BB][^GenPipes_] Tumor-Pair and RNA-seq analytical pipelines (detailed documentation can be found [here][GenPipes_RTD]). Specific parameter values and reference databases are tracked in corresponding *\*.ini* files found under the `parameters` directory. An explanation of the various steps is detailed below.

## WGS
For WGS samples, [GATK][GATK_]’s best practices and procedures were followed. Quality trimmed and adapter-clipped reads were first aligned to the GRCh38 reference[^genome_ref_] with BWA-MEM[^BWA-MEM_]. Alignments were then sorted, realigned around Indels, and marked for duplicates. Base qualities were improved using Base Quality Score Recalibration (BQSR). A MultiQC[^MultiQC_] report was generated per patient to flag any inconsistency in overall coverage, QC bias, tumor purity ([PURPLE][purple_]) and normal or tumor contamination and concordance estimations (ConPair[^ConPair_]).  Somatic and germline calls were generated using an ensemble approach combining four independent variant callers: GATK MuTect2, Strelka2[^Strelka2_], VarDict[^VarDict_] and VarScan2[^VarScan2_]. Somatic and germline variants identified in two or more callers were further annotated and prioritized using the PCGR/CSPR[^PCGR_] reporting system to classify somatic and germline calls using ACMP/AMP classification, perform tumor mutation burden (TMB) estimation, microsatellite instability (MSI) classification, mutational signature estimations, and kataegis detection. Structural Variants were called using GRIDSS with PURPLE and LINX[^GRIDSS_PURPLE_LINX_] as an interpretation tool, and CNAs were called with CNVKit[^CNVKit_]. 

## WTS
For WTS, transcript abundance was estimated using the Kallisto[^Kallisto_] pseudoaligner from quality trimmed and adapter-clipped reads. For WTS variant calling, a full alignment to the same GRCh38 reference was performed using STAR[^STAR_] and the alignments were sorted, realigned and duplicates were marked using GATK best practices and procedures. Fusions in WTS data were assessed using both STAR-FUSION[^STAR-FUSION_] and Arriba[^Arriba_] fusion callers and reported using AnnoFuse[^AnnoFuse_]. Finally, aligned reads were used to call variants with the GATK haplotype caller. Variant calls were then filtered, annotated and prioritized using PCGR/CSPR reporting systems.

In both WGS and WTS, reports were generated with MultiQC for visualization of key quality metrics and allow for manual validation.    

# References
[^GenPipes_]: Bourgey M, Dali R, Eveleigh R, Chen KC, Letourneau L, Fillon J, Michaud M, Caron M, Sandoval J, Lefebvre F, Leveque G, Mercier E, Bujold D, Marquis P, Van PT, Anderson de Lima Morais D, Tremblay J, Shao X, Henrion E, Gonzalez E, Quirion PO, Caron B, Bourque G. GenPipes: an open-source framework for distributed and scalable genomic analyses. Gigascience. 2019 Jun 1;8(6):giz037. doi: 10.1093/gigascience/giz037. PMID: 31185495; PMCID: PMC6559338.

[GenPipes_RTD]: https://genpipes.readthedocs.io


[GenPipes_BB]: https://bitbucket.org/mugqic/genpipes

[GATK_]: https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows

[^genome_ref_]: GCA_000001405.15 no alt analysis set downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz


[^BWA-MEM_]: Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN]


[^MultiQC_]: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PMID: 27312411; PMCID: PMC5039924.


[purple_]: https://github.com/hartwigmedical/hmftools/tree/master/purple

[^ConPair_]: Bergmann EA, Chen BJ, Arora K, Vacic V, Zody MC. Conpair: concordance and contamination estimator for matched tumor-normal pairs. Bioinformatics. 2016 Oct 15;32(20):3196-3198. doi: 10.1093/bioinformatics/btw389. Epub 2016 Jun 26. PMID: 27354699; PMCID: PMC5048070.

[^Strelka2_]: Kim S, Scheffler K, Halpern AL, Bekritsky MA, Noh E, Källberg M, Chen X, Kim Y, Beyter D, Krusche P, Saunders CT. Strelka2: fast and accurate calling of germline and somatic variants. Nat Methods. 2018 Aug;15(8):591-594. doi: 10.1038/s41592-018-0051-x. Epub 2018 Jul 16. PMID: 30013048.

[^VarDict_]: Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, Johnson J, Dougherty B, Barrett JC, Dry JR. VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. Nucleic Acids Res. 2016 Jun 20;44(11):e108. doi: 10.1093/nar/gkw227. Epub 2016 Apr 7. PMID: 27060149; PMCID: PMC4914105.

[^VarScan2_]: Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, Miller CA, Mardis ER, Ding L, Wilson RK. VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Res. 2012 Mar;22(3):568-76. doi: 10.1101/gr.129684.111. Epub 2012 Feb 2. PMID: 22300766; PMCID: PMC3290792.

[^PCGR_]: Nakken S, Fournous G, Vodák D, Aasheim LB, Myklebost O, Hovig E. Personal Cancer Genome Reporter: variant interpretation report for precision oncology. Bioinformatics. 2018 May 15;34(10):1778-1780. doi: 10.1093/bioinformatics/btx817. PMID: 29272339; PMCID: PMC5946881.

[^GRIDSS_PURPLE_LINX_]: GRIDSS, PURPLE, LINX: Unscrambling the tumor genome via integrated analysis of structural variation and copy number. Daniel L. Cameron, Jonathan Baber, Charles Shale, Anthony T. Papenfuss, Jose Espejo Valle-Inclan, Nicolle Besselink, Edwin Cuppen, Peter Priestley. bioRxiv 781013; doi: https://doi.org/10.1101/781013

[^CNVKit_]: Talevich E, Shain AH, Botton T, Bastian BC. CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. PLoS Comput Biol. 2016 Apr 21;12(4):e1004873. doi: 10.1371/journal.pcbi.1004873. PMID: 27100738; PMCID: PMC4839673.

[^Kallisto_]: Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol. 2016 May;34(5):525-7. doi: 10.1038/nbt.3519. Epub 2016 Apr 4. Erratum in: Nat Biotechnol. 2016 Aug 9;34(8):888. PMID: 27043002.

[^STAR_]: Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635. Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.

[^STAR-FUSION_]: STAR-Fusion: Fast and Accurate Fusion Transcript Detection from RNA-Seq. Brian J. Haas, Alex Dobin, Nicolas Stransky, Bo Li, Xiao Yang, Timothy Tickle, Asma Bankapur, Carrie Ganote, Thomas G. Doak, Nathalie Pochet, Jing Sun, Catherine J. Wu, Thomas R. Gingeras, Aviv Regev. bioRxiv 120295; doi: https://doi.org/10.1101/120295

[^Arriba_]: Uhrig S, Ellermann J, Walther T, Burkhardt P, Fröhlich M, Hutter B, Toprak UH, Neumann O, Stenzinger A, Scholl C, Fröhling S, Brors B. Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Res. 2021 Mar;31(3):448-460. doi: 10.1101/gr.257246.119. Epub 2021 Jan 13. PMID: 33441414; PMCID: PMC7919457.

[^AnnoFuse_]: Gaonkar KS, Marini F, Rathi KS, Jain P, Zhu Y, Chimicles NA, Brown MA, Naqvi AS, Zhang B, Storm PB, Maris JM, Raman P, Resnick AC, Strauch K, Taroni JN, Rokita JL. annoFuse: an R Package to annotate, prioritize, and interactively explore putative oncogenic RNA fusions. BMC Bioinformatics. 2020 Dec 14;21(1):577. doi: 10.1186/s12859-020-03922-7. PMID: 33317447; PMCID: PMC7737294."""
    html = markdown.markdown(data, extensions=extensions, extension_configs=extension_configs)
    with open(methods_file, "w", encoding="utf-8", errors="xmlcharrefreplace") as out_file:
        out_file.write(html)

if __name__ == '__main__':
    main()
    #Update db with the objects
