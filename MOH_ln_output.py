#!/usr/bin/env python3

import argparse
import sys
import os
import datetime
import progressbar
import pandas as pd
from  DB_OPS import (
    create_connection,
    extract_sample_metrics,
    extract_sample_details,
    extract_fileloc_field,
    extract_sample_names,
    extract_patient_status
    )

WIDGETS = [' [', progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') - ', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']

def main():
    parser = argparse.ArgumentParser(prog='MOH_ln_output.py', description="Hardlinks files matching criterias into /lustre03/project/6007512/C3G/projects/share/MOH for delivery.")
    parser.add_argument('--black_list', required=False, help="path/to file for patients to be ignored.")
    args = parser.parse_args()

    connection = create_connection("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db")
    # Test String
    # patients = ["MoHQ-CM-3-1-17057-1D","MoHQ-CM-3-2-61285-1D","MoHQ-CM-1-9-2804-1-D"]
    patients = extract_sample_names(connection)

    black_list = []
    if args.black_list:
        with open(args.black_list, "r") as black_list_file:
            for line in black_list_file:
                black_list.append(line.strip())

    patients = list(filter(lambda i: i not in black_list, patients))
    # print(patients)

    # exit()
    with progressbar.ProgressBar(max_value=len(patients), widgets=WIDGETS) as progress:
        for index, patient in enumerate(patients, 1):
            # if patient != "MoHQ-MU-16-4":
            sample = SampleData(connection, patient)
            # print(f"{sample.Sample} {sample.DNA_N} {sample.DNA_T}")
            # Check if samples reach the threashold for delivery
            # Check that DNA_T dedup coverage is over 80
            # Check that DNA_N dedup coverage is over 30
            # Check that processing is complete
            # Check that RNA spots is over 100000000
            # if sample.rna not in ("MoHQ-GC-17-197-OC1-1RT", "MoHQ-GC-17-749-OC1-1RT", "MoHQ-GC-17-997-OC1-1RT"):
            dna = False
            rna = False
            # print(f"\n{sample.dna_n} {sample.dna_t}")
            if  extract_sample_metrics(sample.conn, sample.dna_n, "WGS_Dedup_Coverage") == "NA" or extract_sample_metrics(sample.conn, sample.dna_t, "WGS_Dedup_Coverage") == "NA" or extract_patient_status(sample.conn, sample.sample, "dna_pipeline_execution") == "NA":
                dna = False
            elif float(extract_sample_metrics(sample.conn, sample.dna_n, "WGS_Dedup_Coverage")) > 30 and float(extract_sample_metrics(sample.conn, sample.dna_t, "WGS_Dedup_Coverage")) > 80:
                # print(patient, sample.dna_n, extract_sample_metrics(sample.conn, sample.dna_n, "WGS_Dedup_Coverage"), sample.dna_t, extract_sample_metrics(sample.conn, sample.dna_t, "WGS_Dedup_Coverage"), sep=",")
                dna = True
            if  extract_sample_metrics(sample.conn, sample.rna, "Raw_Reads_Count") == "NA" or extract_patient_status(sample.conn, sample.sample, "rna_pipeline_light_execution") == "NA":
                rna = False
            elif float(extract_sample_metrics(sample.conn, sample.rna, "Raw_Reads_Count")) > 100000000:
                rna = True
                # print(sample.rna)

            # if patient == "MoHQ-JG-9-5":
            #     print(dna)
            # Folders used for Delivery
            # base_folder = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/' # Base Folder
            out_folder = '/lustre03/project/6007512/C3G/projects/share/MOH' # Output Folder
            # Contains Warnings.txt Readme.txt Log.txt and all subfolders
            out_folder = os.path.join(out_folder, sample.institution, sample.cohort, sample.sample_true)
            # Contains raw bams and fastqs
            raw_folder = os.path.join(out_folder, "raw_data")
            # Contains all variants the subfolder
            var_folder = os.path.join(out_folder, "variants")
            # Contains all the vcfs from the callers
            cal_folder = os.path.join(var_folder, "caller_vcfs")
            # Contains raw cnv
            raw_cnv_folder = os.path.join(out_folder, "raw_cnv")
            # Contains the analysis bams
            align_folder = os.path.join(out_folder, "alignment")
            # Contains the ini files
            param_folder = os.path.join(out_folder, "parameters")
            # Contains the reports
            reports_folder = os.path.join(out_folder, "reports")
            pcgr_folder = os.path.join(reports_folder, "pcgr")
            # Contains expression from Kallisto
            expression_folder = os.path.join(out_folder, "expression")

            # updated keeps track of if we need to update the metrics and warnings file
            # old_log stores the data in the log file to see if things need updating
            updated = False
            old_log = {}
            # See if the directory is created and if so check for file updates.
            log = os.path.join(out_folder, "log.txt")
            warn = os.path.join(out_folder, "Warnings.txt")
            if os.path.isdir(out_folder):
                # Load the log file into a dictionary for checking if updates are necessary.
                with open(log, "r") as log_in:
                    #print("logging")
                    for line in log_in:
                        line.rstrip()
                        fields = line.split(",")
                        old_log[fields[0]] = fields[1]
                #log = open(os.path.join(out_folder, "log.txt"), "a")

            elif dna or rna:
                os.makedirs(out_folder)
                os.makedirs(raw_folder)
                os.makedirs(var_folder)
                os.makedirs(raw_cnv_folder)
                os.makedirs(cal_folder)
                os.makedirs(align_folder)
                os.makedirs(param_folder)
                # os.makedirs(tracks_folder)
                os.makedirs(reports_folder)
                os.makedirs(pcgr_folder)
                os.makedirs(expression_folder)

                # Populate the general files
                with open(log, "w") as log_file:
                    log_file.write("File,File_Creation_Date,Date_Added,Details\n")

                generate_readme(out_folder, sample.sample_true, sample.dna_n, sample.dna_t, sample.rna)
                log_new("Readme.txt", log, None, "Created")
            # else:
            #     continue

            # Populate dna data
            if dna:
                beluga_bam_dna_n = extract_fileloc_field(connection, sample.sample, "Beluga_BAM_DNA_N")
                updated = get_link_log(beluga_bam_dna_n, raw_folder, f"{sample.dna_n}.bam", log, updated, old_log)
                beluga_bam_dna_t = extract_fileloc_field(connection, sample.sample, "Beluga_BAM_DNA_T")
                updated = get_link_log(beluga_bam_dna_t, raw_folder, f"{sample.dna_t}.bam", log, updated, old_log)

                dna_vcf_g = extract_fileloc_field(connection, sample.sample, "DNA_VCF_G")
                updated = get_link_log(dna_vcf_g, var_folder, f"{sample.sample_true}.ensemble.germline.vt.annot.vcf.gz", log, updated, old_log)
                dna_vcf_s = extract_fileloc_field(connection, sample.sample, "DNA_VCF_S")
                updated = get_link_log(dna_vcf_s, var_folder, f"{sample.sample_true}.ensemble.somatic.vt.annot.vcf.gz", log, updated, old_log)
                mutect2_germline_vcf = extract_fileloc_field(connection, sample.sample, "Mutect2_Germline_vcf")
                updated = get_link_log(mutect2_germline_vcf, cal_folder, f"{sample.sample_true}.mutect2.germline.vcf.gz", log, updated, old_log)
                mutect2_somatic_vcf = extract_fileloc_field(connection, sample.sample, "Mutect2_Somatic_vcf")
                updated = get_link_log(mutect2_somatic_vcf, cal_folder, f"{sample.sample_true}.mutect2.somatic.vt.vcf.gz", log, updated, old_log)
                strelka2_germline_vcf = extract_fileloc_field(connection, sample.sample, "strelka2_Germline_vcf")
                updated = get_link_log(strelka2_germline_vcf, cal_folder, f"{sample.sample_true}.strelka2.germline.vt.vcf.gz", log, updated, old_log)
                strelka2_somatic_vcf = extract_fileloc_field(connection, sample.sample, "strelka2_Somatic_vcf")
                updated = get_link_log(strelka2_somatic_vcf, cal_folder, f"{sample.sample_true}.strelka2.somatic.vt.vcf.gz", log, updated, old_log)
                vardict_germline_vcf = extract_fileloc_field(connection, sample.sample, "vardict_Germline_vcf")
                updated = get_link_log(vardict_germline_vcf, cal_folder, f"{sample.sample_true}.vardict.germline.vt.vcf.gz", log, updated, old_log)
                vardict_somatic_vcf = extract_fileloc_field(connection, sample.sample, "vardict_Somatic_vcf")
                updated = get_link_log(vardict_somatic_vcf, cal_folder, f"{sample.sample_true}.vardict.somatic.vt.vcf.gz", log, updated, old_log)
                varscan2_germline_vcf = extract_fileloc_field(connection, sample.sample, "varscan2_Germline_vcf")
                updated = get_link_log(varscan2_germline_vcf, cal_folder, f"{sample.sample_true}.varscan2.germline.vt.vcf.gz", log, updated, old_log)
                varscan2_somatic_vcf = extract_fileloc_field(connection, sample.sample, "varscan2_Somatic_vcf")
                updated = get_link_log(varscan2_somatic_vcf, cal_folder, f"{sample.sample_true}.varscan2.somatic.vt.vcf.gz", log, updated, old_log)
                cnvkit_vcf = extract_fileloc_field(connection, sample.sample, "cnvkit_vcf")
                updated = get_link_log(cnvkit_vcf, raw_cnv_folder, f"{sample.sample_true}.cnvkit.vcf.gz", log, updated, old_log)
                # Add md5s and index
                final_dna_bam_n = extract_fileloc_field(connection, sample.sample, "Final_DNA_BAM_N")
                updated = get_link_log(final_dna_bam_n, align_folder, f"{sample.dna_n}.bam", log, updated, old_log)
                if final_dna_bam_n != "NA":
                    final_dna_bam_n_index = final_dna_bam_n + ".bai"
                    if os.path.exists(final_dna_bam_n_index):
                        updated = get_link_log(final_dna_bam_n_index, align_folder, f"{sample.dna_n}.bam.bai", log, updated, old_log)
                    final_dna_bam_n_md5 = final_dna_bam_n + ".md5"
                    if os.path.exists(final_dna_bam_n_md5):
                        updated = get_link_log(final_dna_bam_n_md5, align_folder, f"{sample.dna_n}.bam.md5", log, updated, old_log)
                final_dna_bam_t = extract_fileloc_field(connection, sample.sample, "Final_DNA_BAM_T")
                updated = get_link_log(final_dna_bam_t, align_folder, f"{sample.dna_t}.bam", log, updated, old_log)
                if final_dna_bam_t != "NA":
                    final_dna_bam_t_index = final_dna_bam_t + ".bai"
                    if os.path.exists(final_dna_bam_t_index):
                        updated = get_link_log(final_dna_bam_t_index, align_folder, f"{sample.dna_t}.bam.bai", log, updated, old_log)
                    final_dna_bam_t_md5 = final_dna_bam_t + ".md5"
                    if os.path.exists(final_dna_bam_t_md5):
                        updated = get_link_log(final_dna_bam_t_md5, align_folder, f"{sample.dna_t}.bam.md5", log, updated, old_log)

                dna_multiqc = extract_fileloc_field(connection, sample.sample, "DNA_MultiQC")
                updated = get_link_log(dna_multiqc, reports_folder, f"{sample.sample_true}_D.multiqc.html", log, updated, old_log)
                pcgr_report = extract_fileloc_field(connection, sample.sample, "pcgr_report")
                updated = get_link_log(pcgr_report, reports_folder, f"{sample.sample_true}.pcgr.html", log, updated, old_log)

                pcgr_maf = extract_fileloc_field(connection, sample.sample, "pcgr_maf")
                updated = get_link_log(pcgr_maf, pcgr_folder, f"{sample.sample_true}.acmg.grch38.maf", log, updated, old_log)
                pcgr_snvs_indels = extract_fileloc_field(connection, sample.sample, "pcgr_snvs_indels")
                updated = get_link_log(pcgr_snvs_indels, pcgr_folder, f"{sample.sample_true}.acmg.grch38.snvs_indels.tiers.tsv", log, updated, old_log)
                pcgr_cna_segments = extract_fileloc_field(connection, sample.sample, "pcgr_cna_segments")
                updated = get_link_log(pcgr_cna_segments, pcgr_folder, f"{sample.sample_true}.acmg.grch38.cna_segments.tsv.gz", log, updated, old_log)

                tp_ini = extract_fileloc_field(connection, sample.sample, "TP_ini")
                updated = get_link_log(tp_ini, param_folder, f"{sample.sample_true}.TumourPair.ini", log, updated, old_log)

            if rna:
                # print ("rna_TEST")
                beluga_fastq_1_rna = extract_fileloc_field(connection, sample.sample, "Beluga_fastq_1_RNA")
                updated = get_link_log(beluga_fastq_1_rna, raw_folder, f"{sample.rna}_R1.fastq.gz", log, updated, old_log)
                beluga_fastq_2_rna = extract_fileloc_field(connection, sample.sample, "Beluga_fastq_2_RNA")
                updated = get_link_log(beluga_fastq_2_rna, raw_folder, f"{sample.rna}_R2.fastq.gz", log, updated, old_log)

                rna_abundance = extract_fileloc_field(connection, sample.sample, "RNA_Abundance")
                updated = get_link_log(rna_abundance, expression_folder, f"{sample.rna}.abundance_transcripts.tsv", log, updated, old_log)
                if rna_abundance != "NA":
                    rna_abundance_genes = rna_abundance.replace("transcripts", "genes")
                    if os.path.exists(rna_abundance_genes):
                        updated = get_link_log(rna_abundance_genes, expression_folder, f"{sample.rna}.abundance_genes.tsv", log, updated, old_log)

                # Not implemented yet
                # rna_vcf = extract_fileloc_field(connection, sample.sample, "RNA_VCF")
                # updated = get_link_log(rna_vcf, var_folder, f"{sample.rna}.rna.hc.vcf.gz", log, updated, old_log)

                # Not implemented yet
                # final_rna_bam_variants = extract_fileloc_field(connection, sample.sample, "Final_RNA_BAM_variants")
                # updated = get_link_log(final_rna_bam_variants, align_folder, f"{sample.rna}.variants.bam", log, updated, old_log)
                # if final_rna_bam_variants != "NA":
                #     final_rna_bam_index = final_rna_bam_variants + ".bai"
                #     if os.path.exists(final_rna_bam_index):
                #         updated = get_link_log(final_rna_bam_index, align_folder, f"{sample.rna}.bam.bai", log, updated, old_log)
                #     final_rna_bam_md5 = final_rna_bam_variants + ".md5"
                #     if os.path.exists(final_rna_bam_md5):
                #         updated = get_link_log(final_rna_bam_md5, align_folder, f"{sample.rna}.bam.md5", log, updated, old_log)

                # Not implemented yet
                # rna_multiqc = extract_fileloc_field(connection, sample.sample, "RNA_MultiQC")
                # updated = get_link_log(rna_multiqc, reports_folder, f"{sample.sample_true}_R.multiqc.html", log, updated, old_log)
                # annofuse = extract_fileloc_field(connection, sample.sample, "AnnoFuse")
                # updated = get_link_log(annofuse, reports_folder, f"{sample.rna}.anno_fuse", log, updated, old_log)
                # gridss = extract_fileloc_field(connection, sample.sample, "GRIDSS")
                # updated = get_link_log(gridss, reports_folder, f"{sample.rna}.gridss", log, updated, old_log)

                rna_expression_ini = extract_fileloc_field(connection, sample.sample, "RNA_Abundance_ini")
                updated = get_link_log(rna_expression_ini, param_folder, f"{sample.sample_true}.RNA.Light.ini", log, updated, old_log)
                # Not implemented yet
                # rna_variants_ini = extract_fileloc_field(connection, sample.sample, "RNA_Variants_ini")
                # updated = get_link_log(rna_variants_ini, param_folder, f"{sample.sample_true}.RNA.Variants.ini", log, updated, old_log)

            # Not implemented yet
            # if rna and dna:
            #     final_vcf = extract_fileloc_field(connection, sample.sample, "Final_VCF")
            #     updated = get_link_log(final_vcf, var_folder, f"{sample.sample_true}.vcf.gz", log, updated, old_log)

            # If any updates were made, Delete the old warning file and populate a new one.
            if updated:
                if os.path.exists(warn):
                    os.remove(warn)
                    log_new("Warnings.txt", log, None, "Updated")
                    log_new("Key_metrics.csv", log, None, "Updated")
                else:
                    log_new("Warnings.txt", log, None, "Created")
                    log_new("Key_metrics.csv", log, None, "Updated")

                # Add key metrics table for samples
                metrics = pd.read_sql_query(f'select * from KEY_METRICS where Sample="{sample.dna_n}" or Sample="{sample.dna_t}" or Sample="{sample.rna}"', connection)
                metrics.to_csv(f"{reports_folder}.{sample.sample_true}.Key_metrics.csv", index=False)

                # Add warnings file
                warnings = pd.read_sql_query(f'select Sample,Flags,Fails from KEY_METRICS where Sample="{sample.dna_n}" or Sample="{sample.dna_t}" or Sample="{sample.rna}"', connection)
                warnings.to_csv(warn, index=False)
                with open(warn, 'r+') as warnings_file:
                    content = warnings_file.read()
                    warnings_file.seek(0)
                    warnings_file.write("Below are three columns, \"Flags\" indicates values that may be troublesome while \"Fails\" indicates a point of failure. Data may be useable when marked as \"Flags\", but \"Fails\" marked data should be carefully considered. If nothing is present, this data exceeded all standards.\n Data will not be delivered when coverage is labelled \"Fails\" at this time.\n" + content)
            #log.close()
            progress.update(index)

# Key function. Basiclly updates the logs and links the files based on the database.
def get_link_log(input_file, output_folder, output_file, log, updated, old_log):
    # data = extract_fileloc_field(connection, name, column)
    if input_file != "NA":
        new_time = getime(input_file)
        #Portion for updating files if they have been modified
        if os.path.exists(os.path.join(output_folder, output_file)):
            old_time = old_log[output_file]
            if old_time != new_time:
                os.remove(os.path.join(output_folder, output_file))
                os.link(input_file, os.path.join(output_folder, output_file))
                log_new(output_file, log, new_time, "Updated")
                updated = True
        else:
            os.link(input_file, os.path.join(output_folder, output_file))
            log_new(output_file, log, new_time, "File Added")
            updated = True
    return updated

def getime(path):
    """Finding timestamp from file"""
    date = datetime.datetime.fromtimestamp(os.path.getmtime(path))
    return date.strftime("%Y/%m/%d")

class SampleData:
    def __init__(self, connection, sample):
        data = []
        self.conn = connection
        data = extract_sample_details(connection, sample)
        if len(data) <10:
            raise Exception(f'No database entry for {sample}')
        self.sample = data[1]
        self.sample_true = data[1]
        self.institution = data[2]
        self.cohort = data[3]
        self.dna_n = data[4]
        self.dna_n_true = data[5]
        self.dna_t = data[6]
        self.dna_t_true = data[7]
        self.rna = data[8]
        self.rna_true = data[9]
    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)


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



def generate_readme(location, patient, dna_n, dna_t, rna):
    # Add timestamp
    timestamp = datetime.datetime.today().strftime("%Y/%m/%d")
    data = f"""This directory contains the delivered data for {patient} processed by the Canadian Centre for Computational Genomics.
The data will be updated as it becomes available and as such many files may be missing from RNA or DNA upon initial creation of this directory
Should you have concerns, questions, or suggestions, please contact the analysis team at moh-q@computationalgenomics.ca
Within this directory you will find the results of the analysis for a single patient contained in 7 subdirectories and three files:

Readme.txt      This file
log.txt         Log file containing the dates of transfers and if files have been updated
Warnings.txt    Contains details of any warnings and whether they caused a failure of this analysis

raw_data/               Contains all of the bam's/fastqs from the sequencer. BAM files here include both mapped and unmapped reads and can be converted to the FASTQ format with tools such as SamToFastq.
    {dna_n}.bam         Raw DNA reads for the Normal sample
    {dna_t}.bam         Raw DNA reads for the Tumor sample
    {rna}.fastq         Raw RNA reads for the Tumor sample

variants/                                       Contains the vcfs related to variant calls
    {patient}.ensemble.germline.vt.annot.vcf.gz         Germline Variants found in any of the callers
    {patient}.ensemble.somatic.vt.annot.vcf.gz          Somatic Variants found in any of the callers
*   {patient}.rna.hc.vcf.gz                             Variants found using RNA sample (/!\\ Not yet available)
*   {patient}.vcf.gz                                    Contains the results of all callers for both DNA and RNA (/!\\ Not yet available)
    caller_vcfs/                                Contains the vcfs produced from individual callers on the DNA samples
        {patient}.mutect2.germline.vcf.gz              Germline results for mutect2
        {patient}.mutect2.somatic.vt.vcf.gz            Somatic results for mutect2
        {patient}.strelka2.germline.vt.vcf.gz          Germline results for strelka2
        {patient}.strelka2.somatic.vt.vcf.gz           Somatic results for strelka2
        {patient}.vardict.germline.vt.vcf.gz           Germline results for vardict
        {patient}.vardict.somatic.vt.vcf.gz            Somatic results for vardict
        {patient}.varscan2.germline.vt.vcf.gz          Germline results for varscan2
        {patient}.varscan2.somatic.vt.vcf.gz           Somatic results for varscan2

raw_cnv/                Contains the raw copy number calls for each patient DNA
    {patient}.cnvkit.vcf.gz     Raw cnvkit output

alignment/              Contains the alignment data for each sample
    {dna_n}.bam                 Alignment of normal against the reference
    {dna_n}.bam.bai             Index of Alignment of normal against the reference
    {dna_t}.bam                 Alignment of tumor against the reference
    {dna_t}.bam.bai             Index of Alignment of tumor against the reference
*   {rna}.variants.bam          Alignment of tumor RNA against the reference used in variants analysis (/!\\ Not yet available)
*   {rna}.variants.bam.bai      Index of Alignment of tumor RNA against the reference used in variants analysis (/!\\ Not yet available)

expression/             Contains the transcripts abundance estimation from Kallisto
    {rna}.abundance_transcripts.tsv       Table with transcript abundance feom Kallisto

reports/                Contains the reports for the experiment
    {patient}_D.multiqc.html    QC report for the DNA analysis
*   {patient}_R.multiqc.html    QC report for the RNA analysis (/!\\ Not yet available)
    {patient}.pcgr.html         Personal Cancer Genome Reporter report; Cf. https://pcgr.readthedocs.io/en/latest
*   {rna}.anno_fuse             Report for fusions detected using RNA (/!\\ Not yet available)
*   {patient}.gridss            Annotated structural variant calls with GRIDSS (/!\\ Not yet available)
*   {patient}.Key_metrics.csv   Metrics used to determine whether the analysis was successful (/!\\ Not yet available)
    pcgr/               Contains raw tables used to generate PCGR report
        {patient}.acmg.grch38.maf
        {patient}.acmg.grch38.snvs_indels.tiers.tsv
        {patient}.acmg.grch38.cna_segments.tsv.gz

parameters/             Contains the records of all the Parameters used in the pipeline analysis
    tumorPair.config.trace.ini              Parameters used in the tumor pair analysis
    RNAseq.expression.config.trace.ini      Parameters used in the RNA expression analysis
*   RNAseq.variants.config.trace.ini        Parameters used in the RNA variant analysis (/!\\ Not yet available)

Generated {timestamp}"""
    with open(os.path.join(location, "Readme.txt"), "w") as readme_file:
        readme_file.write(data)

if __name__ == '__main__':
    main()
    #Update db with the objects
