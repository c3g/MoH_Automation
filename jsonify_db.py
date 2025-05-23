#!/usr/bin/env python3

import json
import glob
import re
import csv
import os
import hashlib
import logging
from datetime import datetime
import numpy as np
import h5py

logging.basicConfig(format='%(levelname)s: %(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    """
    Main function to run the script.
    """
    # prefix_path = "/Users/pstretenowich/Mount_points/beluga"
    prefix_path = "/lustre03/project/6007512"
    main_raw_reads_folder = os.path.join(prefix_path, "C3G/projects/MOH_PROCESSING/MAIN/raw_reads")
    # patient_dict = get_patient_dict("CheckBeluga.txt") #For the recraft of only the missing runs 08/2024 . Adapt acordg to the input
    patient_dict = get_patient_dict(main_raw_reads_folder)
    _, sample_dict = jsonify_run_processing(patient_dict, prefix_path)
    jsonify_transfer(sample_dict, prefix_path)
    jsonify_genpipes_tumourpair_ensemble(sample_dict, prefix_path)
    jsonify_genpipes_rnaseqlight(sample_dict, prefix_path)
    jsonify_genpipes_tumourpair_ensemble(sample_dict,prefix_path)
    jsonify_genpipes_tumourpair_sv(sample_dict, prefix_path)
    jsonify_genpipes_rnaseqlight(sample_dict, prefix_path)
    jsonify_genpipes_rnaseq_cancer(sample_dict, prefix_path)


def jsonify_run_processing(patient_dict, prefix_path):
    """
    Create a json file for each run in the run_metrics folder.
    Return: dict, dict
    """
    readset_dict = {}
    sample_dict = {}
    run_metrics_csv = os.path.join(prefix_path, "C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/*.csv")
    # run_metrics_csv = "/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/220128_A00266_0620_BH527LDSX3_MOHRun01_Q008401-novaseq-run.align_bwa_mem.csv"
    for run_file in glob.glob(run_metrics_csv):
        run_list = []
        with open(run_file, 'rt') as run_file_in:
            reader = csv.DictReader(run_file_in)
            for row in reader:
                run_list.append(row)
                json_output = {
                "operation_platform": "abacus",
                "project_fms_id": None,
                "project_name": "MOH-Q",
                "run_fms_id": None,
                "run_name": f"{run_list[0]['Processing Folder Name']}",
                "run_instrument": "novaseq",
                "run_date": f"{datetime.strptime(run_list[0]['Processing Folder Name'][0:6], '%y%m%d')}",
                "patient": []
                }
        for run_row in run_list:
            sample = run_row['Sample Name']
            if sample.startswith("MoHQ"):
                result = re.search(r"^((MoHQ-(JG|CM|GC|IQ|HM|MU|MR|XX|CQ)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
                patient = result.group(1)
                try:
                    existing_patient = None
                    for position, element in enumerate(json_output["patient"]):
                        if f"{patient_dict[patient][sample][0]}" == element["patient_name"]:
                            existing_patient = element
                            existing_patient_pos = position
                    if existing_patient:
                        patient_json = existing_patient
                        del json_output["patient"][existing_patient_pos]
                    else:
                        patient_json = {
                            "patient_fms_id": None,
                            "patient_name": f"{patient_dict[patient][sample][0]}",
                            "patient_cohort": f"{patient_dict[patient][sample][3]}",
                            "patient_institution": f"{patient_dict[patient][sample][4]}",
                            "sample": []
                            }

                    existing_sample = None
                    sample_name = f"{patient_dict[patient][sample][1]}"
                    for _, element in enumerate(patient_json["sample"]):
                        if sample_name == element["sample_name"]:
                            existing_sample = element
                    if existing_sample:
                        sample_json = existing_sample
                    else:
                        sample_tumour = patient_dict[patient][sample][2].endswith("T")
                        sample_json = {
                            "sample_fms_id": None,
                            "sample_name": sample_name,
                            "sample_tumour": sample_tumour,
                            "readset": []
                            }

                    transfer_folder = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*')
                    # transfer_folder = 'transfer/*'
                    fastq1 = fastq2 = ""
                    for filename in glob.glob(transfer_folder):
                        with open(filename, 'r') as file:
                            for line in file:
                                if patient_dict[patient][sample][1] in line and run_row['Processing Folder Name'] in line:
                                    fields = line.split(",")
                                    if ".bam" in line:
                                        file_json = [
                                            {
                                                "location_uri": f"abacus://{fields[0]}",
                                                "file_name": f"{os.path.basename(fields[0])}",
                                                "file_deliverable": True
                                                }
                                            ]
                                        break
                                    if ".fastq" in line:
                                        if "R1.fastq" in line:
                                            fastq1 = os.path.basename(fields[0])
                                            fastq1_location_uri = fields[0]
                                        elif "R2.fastq" in line:
                                            fastq2 = os.path.basename(fields[0])
                                            fastq2_location_uri = fields[0]
                                        if fastq1 and fastq2:
                                            file_json = [
                                            {
                                                "location_uri": f"abacus://{fastq1_location_uri}",
                                                "file_name": f"{fastq1}",
                                                "file_extra_metadata": {"read_type": "R1"},
                                                "file_deliverable": True
                                                },
                                            {
                                                "location_uri": f"abacus://{fastq2_location_uri}",
                                                "file_name": f"{fastq2}",
                                                "file_extra_metadata": {"read_type": "R2"},
                                                "file_deliverable": True
                                                }
                                            ]
                                            break
                    if not run_row['Clusters']:
                        raw_reads_count_flag = "MISSING"
                    if run_row['Clusters'] =='0':
                        raw_reads_count_flag = "FAILED"
                    else:
                        raw_reads_count_flag = "PASS"
                    if run_row['Library Type'] == "RNASeq":
                        raw_reads_count_flag = rna_raw_reads_count_check(run_row['Clusters'])
                    raw_duplication_rate_flag = "PASS"
                    if run_row['Library Type'] != "RNASeq":
                        raw_duplication_rate_flag = dna_raw_duplication_rate_check(run_row['Dup. Rate (%)'])
                    raw_median_insert_size_flag = median_insert_size_check(run_row['Mapped Insert Size (median)'])
                    raw_mean_insert_size_flag = "PASS"
                    raw_mean_coverage_flag = "PASS"
                    if run_row['Library Type'] != "RNASeq":
                        raw_mean_coverage_flag = dna_raw_mean_coverage_check(run_row['Mean Coverage'], sample_tumour)
                    metric_json = [
                        {
                            "metric_name": "raw_reads_count",
                            "metric_value": f"{run_row['Clusters']}",
                            "metric_flag": raw_reads_count_flag,
                            "metric_deliverable": True
                            },
                        {
                            "metric_name": "raw_duplication_rate",
                            "metric_value": f"{run_row['Dup. Rate (%)']}",
                            "metric_flag": raw_duplication_rate_flag
                            },
                        {
                            "metric_name": "raw_median_insert_size",
                            "metric_value": f"{run_row['Mapped Insert Size (median)']}",
                            "metric_flag": raw_median_insert_size_flag
                            },
                        {
                            "metric_name": "raw_mean_insert_size",
                            "metric_value": f"{run_row['Mapped Insert Size (mean)']}",
                            "metric_flag": raw_mean_insert_size_flag
                            },
                        {
                            "metric_name": "raw_mean_coverage",
                            "metric_value": f"{run_row['Mean Coverage']}",
                            "metric_flag": raw_mean_coverage_flag
                            }
                        ]
                    readset_name = f"{run_row['Sample Name']}.{run_row['Run ID']}_{run_row['Lane']}"
                    readset_dict[readset_name] = (patient, sample)
                    readset_json = {
                        "experiment_sequencing_technology": None,
                        "experiment_type": f"{run_row['Library Type']}",
                        "experiment_library_kit": None,
                        "experiment_kit_expiration_date": None,
                        "readset_name": readset_name,
                        "readset_lane": f"{run_row['Lane']}",
                        "readset_adapter1": f"{run_row['i7 Adapter Sequence']}",
                        "readset_adapter2": f"{run_row['i5 Adapter Sequence']}",
                        "readset_sequencing_type": f"{run_row['Run Type']}",
                        "readset_quality_offset": "33",
                        "file": file_json,
                        "metric": metric_json
                        }
                    sample_json["readset"].append(readset_json)
                    if sample_name in sample_dict:
                        sample_dict[sample_name].append((patient, readset_name))
                    else:
                        sample_dict[sample_name] = [(patient, readset_name)]
                    patient_json["sample"].append(sample_json)
                    json_output["patient"].append(patient_json)
                    # for sample in patient_json["sample"]:
                    #     print(sample, "\n")
                    #     print(sample["sample_name"], "\n")
                except KeyError:
                    pass
        with open(f"jsons/{run_list[0]['Processing Folder Name']}.json", 'w', encoding='utf-8') as f:
            json.dump(json_output, f, ensure_ascii=False, indent=4)

    return readset_dict, sample_dict

def jsonify_transfer(sample_dict, prefix_path):
    """
    Create a json file for each transfer log file in the transfer folder.
    """
    transfer_folder = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*')
    for filename in glob.glob(transfer_folder):
        transfer_dict = {}
        json_output = {
            "operation_platform": "beluga",
            "operation_cmd_line": f"globus transfer --submission-id $sub_id --label $label --batch /lb/project/mugqic/projects/MOH/TEMP/{os.path.basename(filename)} 6c66d53d-a79d-11e8-96fa-0a6d4e044368 278b9bfe-24da-11e9-9fa2-0a06afd4a22e",
            "readset": []
            }
        with open(filename, 'r') as file:
            for line in file:
                if line.startswith("/"):
                    fields = line.split(",")
                    sample = os.path.dirname(fields[1])
                    if sample in sample_dict:
                        abacus_uri = "abacus://"
                        beluga_uri = "beluga:///lustre03/project/6007512/C3G/projects/MOH_PROCESSING/raw_reads/" + sample
                        if ".bam" in line:
                            readset_name = sample + "." + fields[0].split("/")[11].replace("run", "")
                            bam_src_location_uri = abacus_uri + fields[0]
                            bam_dest_location_uri = os.path.join(beluga_uri, os.path.basename(fields[1].strip()))
                            if readset_name in transfer_dict:
                                transfer_dict[readset_name]["bam_dest_location_uri"] = bam_dest_location_uri
                            else:
                                transfer_dict[readset_name] = {
                                    "bam_dest_location_uri": bam_dest_location_uri
                                }
                            transfer_dict[readset_name]["bam_src_location_uri"] = bam_src_location_uri

                        elif ".fastq" in line:
                            run_name = fields[0].split("/")[7]
                            readset_name = sample + "." + run_name.split("_")[1] + "_" + run_name.split("_")[2] + "_" + fields[0].split("/")[8].split(".")[1]
                            if "R1.fastq" in line:
                                fastq1_src_location_uri = abacus_uri + fields[0]
                                fastq1_dest_location_uri = os.path.join(beluga_uri, os.path.basename(fields[1].strip()))
                                if readset_name in transfer_dict:
                                    transfer_dict[readset_name]["fastq1_dest_location_uri"] = fastq1_dest_location_uri
                                else:
                                    transfer_dict[readset_name] = {
                                        "fastq1_dest_location_uri": fastq1_dest_location_uri
                                    }
                                transfer_dict[readset_name]["fastq1_src_location_uri"] = fastq1_src_location_uri
                            elif "R2.fastq" in line:
                                fastq2_src_location_uri = abacus_uri + fields[0]
                                fastq2_dest_location_uri = os.path.join(beluga_uri, os.path.basename(fields[1].strip()))
                                if readset_name in transfer_dict:
                                    transfer_dict[readset_name]["fastq2_dest_location_uri"] = fastq2_dest_location_uri
                                else:
                                    transfer_dict[readset_name] = {
                                        "fastq2_dest_location_uri": fastq2_dest_location_uri
                                    }
                                transfer_dict[readset_name]["fastq2_src_location_uri"] = fastq2_src_location_uri

        for readset in transfer_dict:
            if "bam_src_location_uri" in transfer_dict[readset]:
                file_json = [
                    {
                        "src_location_uri": transfer_dict[readset]["bam_src_location_uri"],
                        "dest_location_uri": transfer_dict[readset]["bam_dest_location_uri"]
                    }
                ]
            elif "fastq1_src_location_uri" in transfer_dict[readset]:
                file_json = [
                    {
                        "src_location_uri": transfer_dict[readset]["fastq1_src_location_uri"],
                        "dest_location_uri": transfer_dict[readset]["fastq1_dest_location_uri"]

                    },
                    {
                        "src_location_uri": transfer_dict[readset]["fastq2_src_location_uri"],
                        "dest_location_uri": transfer_dict[readset]["fastq2_dest_location_uri"],

                    }
                ]
            else:
                file_json = []
            readset_json = {
                "readset_name": readset,
                "file": file_json,
                }
            json_output["readset"].append(readset_json)
            # print(json.dumps(readset_json, indent=4))
        if transfer_dict:
            with open(f"jsons/{os.path.basename(filename).replace('.txt', '.json')}", 'w', encoding='utf-8') as f:
                json.dump(json_output, f, ensure_ascii=False, indent=4)


def jsonify_genpipes_tumourpair_sv(sample_dict, prefix_path):
    """
    Create a json file for each patient in the sample_dict for the TumorPair SV pipeline.
    """
    ini_file = os.path.join(prefix_path, "C3G/projects/MOH_PROCESSING/MAIN/TumorPair.config.trace.ini")
    to_parse = False
    ini_content = []
    with open(ini_file, 'r') as file:
        for line in file:
            if "base.ini" in line:
                genpipes_version = re.findall(r"/genpipes-.+?/", line)[0].split("-")[-1][:-1]
            else:
                genpipes_version = "unknown"
            if "[DEFAULT]" in line:
                to_parse = True
            if to_parse:
                ini_content.append(line)
    # From where to parse genpipes version and cmd line
    json_output = {
        "project_name": "MOH-Q",
        "operation_config_name": "genpipes_ini",
        "operation_config_version": f"{genpipes_version}",
        "operation_config_md5sum": f"{md5(ini_file)}",
        "operation_config_data": f"{''.join(ini_content)}",
        "operation_platform": "beluga",
        "operation_cmd_line": "tumor_pair.py -s 1-13,15-38 -j slurm -r readset.txt -c $MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.base.ini $MUGQIC_PIPELINES_HOME/pipelines/common_ini/beluga.ini $MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.extras.ini $MUGQIC_PIPELINES_HOME/resources/genomes/config/Homo_sapiens.GRCh38.ini Custom_ini/Custom_MOH_topups.ini -g TP_run.sh",
        "operation_name": "GenPipes_TumorPair.sv",
        "sample": []
        }
    sample_dict_dna = {}
    for sample in sample_dict:
        if not sample.endswith("RT"):
            sample_dict_dna[sample] = sample_dict[sample]
    length = len(sample_dict_dna)
    i = 0
    for sample in sample_dict_dna:
        i += 1
        eta = round(float(100*i/length), 2)
        logger.info(f"jsonify_genpipes_tumourpair_sv - {eta}%")
        patient = sample_dict_dna[sample][0][0]
        tumour = sample.endswith("DT")
        sample_json = {
            "sample_name": f"{sample}",
            "readset": []
        }
        job_jsons = []
        job_jsons.append(gridss_paired_somatic(patient, sample, prefix_path))
        job_jsons.append(gripss_filter_somatic(patient, sample, prefix_path))
        job_jsons.append(gripss_filter_germline(patient, sample, prefix_path))
        job_jsons.append(purple_purity_sv(patient, sample, prefix_path))

        if tumour:
            job_jsons.append(linx_annotations_germline(patient, sample, prefix_path))
            job_jsons.append(linx_annotations_somatic(patient, sample, prefix_path))
            job_jsons.append(linx_plot(patient, prefix_path))

        for (_, readset) in sample_dict_dna[sample]:
            readset_json = {
                "readset_name": f"{readset}",
                "job": []
            }
            job_jsons = list(filter(lambda job_json: job_json is not None, job_jsons))
            for job_json in job_jsons:
                readset_json["job"].append(job_json)
            if readset_json["job"]:
                sample_json["readset"].append(readset_json)
        if sample_json["readset"]:
            json_output["sample"].append(sample_json)

    if json_output["sample"]:
        with open("jsons/genpipes_tumourpair_sv.json", 'w', encoding='utf-8') as f:
            json.dump(json_output, f, ensure_ascii=False, indent=4)

def jsonify_genpipes_tumourpair_ensemble(sample_dict, prefix_path):
    """
    Create a json file for each patient in the sample_dict for the TumorPair ensemble pipeline.
    """
    ini_file = os.path.join(prefix_path, "C3G/projects/MOH_PROCESSING/MAIN/TumorPair.config.trace.ini")
    to_parse = False
    ini_content = []
    with open(ini_file, 'r') as file:
        for line in file:
            if "base.ini" in line:
                genpipes_version = re.findall(r"/genpipes-.+?/", line)[0].split("-")[-1][:-1]
            else:
                genpipes_version = "unknown"
            if "[DEFAULT]" in line:
                to_parse = True
            if to_parse:
                ini_content.append(line)
    # From where to parse genpipes version and cmd line
    json_output = {
        "project_name": "MOH-Q",
        "operation_config_name": "genpipes_ini",
        "operation_config_version": f"{genpipes_version}",
        "operation_config_md5sum": f"{md5(ini_file)}",
        "operation_config_data": f"{''.join(ini_content)}",
        "operation_platform": "beluga",
        "operation_cmd_line": "tumor_pair.py -s 1-13,15-38 -j slurm -r readset.txt -c $MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.base.ini $MUGQIC_PIPELINES_HOME/pipelines/common_ini/beluga.ini $MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.extras.ini $MUGQIC_PIPELINES_HOME/resources/genomes/config/Homo_sapiens.GRCh38.ini Custom_ini/Custom_MOH_topups.ini -g TP_run.sh",
        "operation_name": "GenPipes_TumorPair.ensemble",
        "sample": []
        }
    sample_dict_dna = {}
    for sample in sample_dict:
        if not sample.endswith("RT"):
            sample_dict_dna[sample] = sample_dict[sample]
    length = len(sample_dict_dna)
    i = 0
    for sample in sample_dict_dna:
        i += 1
        eta = round(float(100*i/length), 2)
        logger.info(f"jsonify_genpipes_tumourpair - {eta}%")
        patient = sample_dict_dna[sample][0][0]
        tumour = sample.endswith("DT")
        sample_json = {
            "sample_name": f"{sample}",
            "readset": []
        }
        job_jsons = []

        job_jsons.append(sym_link_final_bam_tumourpair(patient, sample, tumour, prefix_path))
        job_jsons.append(gatk_variant_annotator_germline_tumourpair(patient, prefix_path))
        job_jsons.append(gatk_variant_annotator_somatic_tumourpair(patient, prefix_path))
        job_jsons.append(paired_mutect2_tumourpair(patient, prefix_path))
        job_jsons.append(vardict_paired_tumourpair(patient, prefix_path))
        job_jsons.append(paired_varscan2_tumourpair(patient, prefix_path))
        job_jsons.append(cnvkit_batch_tumourpair(patient, prefix_path))
        job_jsons.append(recalibration_tumourpair(patient,sample, prefix_path))
        job_jsons.append(strelka2_paired_germline_tumourpair(patient, prefix_path))
        job_jsons.append(report_pcgr_tumourpair(patient, prefix_path))
        job_jsons.append(tumour_pair_multiqc(patient, prefix_path))

        # Conpair
        job_jsons.append(extract_conpair(patient, sample, tumour, prefix_path))
        # Purple
        job_jsons.append(extract_purple(sample, patient, prefix_path))
        # Qualimap
        job_jsons.append(extract_qualimap(patient,sample, prefix_path))
        # Picard
        job_jsons.append(extract_picard_tumourpair(patient,sample, prefix_path))
        for (_, readset) in sample_dict_dna[sample]:
            readset_json = {
                "readset_name": f"{readset}",
                "job": []
            }
            job_jsons = list(filter(lambda job_json: job_json is not None, job_jsons))
            for job_json in job_jsons:
                readset_json["job"].append(job_json)
            if readset_json["job"]:
                sample_json["readset"].append(readset_json)
        if sample_json["readset"]:
            json_output["sample"].append(sample_json)
    if json_output["sample"]:
        with open("jsons/genpipes_tumourpair_en.json", 'w', encoding='utf-8') as f:
            json.dump(json_output, f, ensure_ascii=False, indent=4)


def jsonify_genpipes_rnaseqlight(sample_dict, prefix_path):
    """
    Create a json file for each patient in the sample_dict for the RnaSeqLight pipeline.
    """
    ini_file = os.path.join(prefix_path, "C3G/projects/MOH_PROCESSING/MAIN/RnaSeqLight.config.trace.ini")
    to_parse = False
    ini_content = []
    with open(ini_file, 'r') as file:
        for line in file:
            if "base.ini" in line:
                genpipes_version = re.findall(r"/genpipes-.+?/", line)[0].split("-")[-1][:-1]
            else:
                genpipes_version = "unknown"
            if "[DEFAULT]" in line:
                to_parse = True
            if to_parse:
                ini_content.append(line)
    # From where to parse genpipes version and cmd line
    json_output = {
        "project_name": "MOH-Q",
        "operation_config_name": "genpipes_ini",
        "operation_config_version": f"{genpipes_version}",
        "operation_config_md5sum": f"{md5(ini_file)}",
        "operation_config_data": f"{''.join(ini_content)}",
        "operation_platform": "beluga",
        "operation_cmd_line": "rnaseq_light.py -s 1-4 -j slurm -r readset.txt -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq_light/rnaseq_light.base.ini $MUGQIC_PIPELINES_HOME/pipelines/common_ini/beluga.ini $MUGQIC_PIPELINES_HOME/resources/genomes/config/Homo_sapiens.GRCh38.ini RNA_light.custom.ini -g RNASeq_light_run.sh",
        "operation_name": "GenPipes_RnaSeqLight",
        "sample": []
        }
    sample_dict_rna = {}
    for sample in sample_dict:
        if sample.endswith("RT"):
            sample_dict_rna[sample] = sample_dict[sample]
    for sample in sample_dict_rna:
        sample_json = {
            "sample_name": f"{sample}",
            "readset": []
        }
        job_jsons = []

        # Picard
        job_jsons.append(extract_picard_rna(sample, prefix_path))
        # Kallisto
        job_jsons.append(extract_kallisto_rnaseqlight(sample, prefix_path))
        for (_, readset) in sample_dict_rna[sample]:
            readset_json = {
                "readset_name": f"{readset}",
                "job": []
            }
            job_jsons = list(filter(lambda job_json: job_json is not None, job_jsons))
            for job_json in job_jsons:
                readset_json["job"].append(job_json)
            if readset_json["job"]:
                sample_json["readset"].append(readset_json)
        if sample_json["readset"]:
            json_output["sample"].append(sample_json)
    if json_output["sample"]:
        with open("jsons/genpipes_rnaseqlight.json", 'w', encoding='utf-8') as f:
            json.dump(json_output, f, ensure_ascii=False, indent=4)

def jsonify_genpipes_rnaseq_cancer(sample_dict, prefix_path):
    """
    Create a json file for each patient in the sample_dict for the RnaSeqCancer pipeline.
    """
    ini_file = os.path.join(prefix_path, "C3G/projects/MOH_PROCESSING/MAIN/RnaSeq.config.trace.ini")
    to_parse = False
    ini_content = []
    with open(ini_file, 'r') as file:
        for line in file:
            if "base.ini" in line:
                genpipes_version = re.findall(r"/genpipes-.+?/", line)[0].split("-")[-1][:-1]
            else:
                genpipes_version = "unknown"
            if "[DEFAULT]" in line:
                to_parse = True
            if to_parse:
                ini_content.append(line)
    # From where to parse genpipes version and cmd line
    json_output = {
        "project_name": "MOH-Q",
        "operation_config_name": "genpipes_ini",
        "operation_config_version": f"{genpipes_version}",
        "operation_config_md5sum": f"{md5(ini_file)}",
        "operation_config_data": f"{''.join(ini_content)}",
        "operation_platform": "beluga",
        "operation_cmd_line": "rnaseq.py -s 1-7,10-23 -j slurm -r readset.txt -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini $MUGQIC_PIPELINES_HOME/pipelines/common_ini/beluga.ini /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/Custom_ini/tumor_rna.moh.ini $MUGQIC_PIPELINES_HOME/resources/genomes/config/Homo_sapiens.GRCh38.ini -g RNASeq_run.sh -t cancer",
        "operation_name": "GenPipes_RnaSeq.cancer",
        "sample": []
        }
    sample_dict_rna = {}
    for sample in sample_dict:
        if sample.endswith("RT"):
            sample_dict_rna[sample] = sample_dict[sample]
        # if sample == "MoHQ-MU-17-1251-OC1-1DT":
        #     sample_dict_rna[sample] = sample_dict[sample]
    for sample in sample_dict_rna:
        sample_json = {
            "sample_name": f"{sample}",
            "readset": []
        }
        job_jsons = []

        # Rnaseqc2
        job_jsons.append(extract_rnaseqc2_rna(sample, prefix_path))
        job_jsons.append(filter_gatk_rnaseqc2(sample, prefix_path))
        job_jsons.append(recalibration_rnaseqc2(sample, prefix_path))
        job_jsons.append(rnaseqc2_multiqc(sample, prefix_path))
        job_jsons.append(report_pcgr_rnaseqc2(sample, prefix_path))
        for (_, readset) in sample_dict_rna[sample]:
            readset_json = {
                "readset_name": f"{readset}",
                "job": []
            }
            job_jsons = list(filter(lambda job_json: job_json is not None, job_jsons))
            for job_json in job_jsons:
                readset_json["job"].append(job_json)
            if readset_json["job"]:
                sample_json["readset"].append(readset_json)
        if sample_json["readset"]:
            json_output["sample"].append(sample_json)
    if json_output["sample"]:
        with open("jsons/genpipes_rnaseqcancer.json", 'w', encoding='utf-8') as f:
            json.dump(json_output, f, ensure_ascii=False, indent=4)

def md5(fname):
    """
    Calculate the md5 hash of a file.
    Return: str
    """
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def get_patient_dict(main_raw_reads_folder):
    """
    Create a dictionary with patient as key and a dictionary with sample as key and a tuple with patient, sample, sample type, cohort and institution as value.
    Return: dict
    """
    patient_dict = {}
    for sample in os.listdir(main_raw_reads_folder):
        if sample.startswith("MoHQ"):
            result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX|CQ)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
            patient = result.group(1)
            cohort = result.group(2)
            institution = result.group(3)
            sample_type = sample[-2:]
            try:
                patient_dict[patient][sample] = (patient, sample, sample_type, cohort, institution)
            except KeyError:
                patient_dict[patient] = {sample: (patient, sample, sample_type, cohort, institution)}
    return patient_dict

def extract_conpair(patient, sample, tumour, prefix_path):
    """
    Extract conpair job, files and metrics from the job_output/conpair_concordance_contamination/conpair_concordance_contamination.pileup.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/conpair_concordance_contamination/conpair_concordance_contamination.pileup.*{sample}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/conpair_concordance_contamination/conpair_concordance_contamination.pileup.*{sample}*.o file found")
    job_json_conpair = {
        "job_name": f"conpair.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": [],
        "metric": []
    }
    # Contamination
    no_contamination = False
    # The file is named after tumour sample only
    contamination_file = glob.glob(os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/metrics', patient + '-*DT.contamination.tsv'))
    # Test if unsure about finding more than 1 contamination file for normal sample based on the glob above.
    # It has to print nothing to be ok, otherwise it means a manual check is required.
    if len(contamination_file) > 1:
        print(f" WARNING: Manual check reauired for patient {patient} as more than 1 contamination file is found: {contamination_file}")
    try:
        with open(contamination_file[0], 'r', encoding="utf-8") as file:
            job_json_conpair["job_status"] = "COMPLETED"
            for line in file:
                if line.startswith('Normal') and not tumour:
                    value = line.split(" ")[-1][:-2]
                    if float(value)>5:
                        flag = "FAILED"
                    else:
                        flag = "PASS"
                    job_json_conpair["metric"].append({
                        "metric_name": "contamination",
                        "metric_value": f"{value}",
                        "metric_flag": f"{flag}",
                        "metric_deliverable": True
                        })
                    job_json_conpair["file"].append({
                        "location_uri": f"beluga://{contamination_file[0]}",
                        "file_name": f"{os.path.basename(contamination_file[0])}"
                        })
                elif line.startswith('Tumor') and tumour:
                    value = line.split(" ")[-1][:-2]
                    if float(value)>5:
                        flag = "FAILED"
                    else:
                        flag = "PASS"
                    job_json_conpair["metric"].append({
                        "metric_name": "contamination",
                        "metric_value": f"{value}",
                        "metric_flag": f"{flag}",
                        "metric_deliverable": True
                        })
                    job_json_conpair["file"].append({
                        "location_uri": f"beluga://{contamination_file[0]}",
                        "file_name": f"{os.path.basename(contamination_file[0])}"
                        })
    except (FileNotFoundError, IndexError):
        no_contamination = True
    # Concordance
    no_concordance = False
    # The file is named after tumour sample only
    filename = glob.glob(os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/metrics', f"{patient}-*DT.concordance.tsv"))
    # Test if unsure about finding more than 1 concordance file for normal sample based on the glob above.
    # It has to print nothing to be ok, otherwise it means a manual check is required.
    if len(filename) > 1:
        print(f" WARNING: Manual check reauired for patient {patient} as more than 1 concordance file is found: {filename}")
    try:
        with open(filename[0], 'r', encoding="utf-8") as file:
            job_json_conpair["job_status"] = "COMPLETED"
            for line in file:
                if line.startswith('Concordance'):
                    value = line.split(" ")[-1][:-2]
                    if float(value)<99:
                        flag = "FAILED"
                    else:
                        flag = "PASS"
                    job_json_conpair["metric"].append({
                        "metric_name": "concordance",
                        "metric_value": f"{value}",
                        "metric_flag": f"{flag}",
                        "metric_deliverable": True
                        })
                    job_json_conpair["file"].append({
                        "location_uri": f"beluga://{filename[0]}",
                        "file_name": f"{os.path.basename(filename[0])}"
                        })
    except (FileNotFoundError, IndexError):
        no_concordance = True
    if no_concordance and no_contamination:
        job_json_conpair = None
    return job_json_conpair


def extract_purple(sample, patient, prefix_path):
    """
    Extract purple job, files and metrics from the job_output/purple/purple.purity.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/purple/purple.purity.*{patient}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/purple/purple.purity.*{patient}*.o file found")
    job_json = {
        "job_name": f"purple.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": [],
    }
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, f'{patient}.strelka2.somatic.purple.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, f'{patient}.strelka2.somatic.purple.vcf.gz.tbi')
    job_json = add_output_file(filename, job_json)
    try:
        filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, 'purple', sample + '.purple.purity.tsv')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            line = lines[1]
            fields = line.split("\t")
            value = float(fields[0])*100
            if float(value)<30:
                flag = "FAILED"
            else:
                flag = "PASS"
            try:
                job_json["metric"].append({
                    "metric_name": "purity",
                    "metric_value": f"{value}",
                    "metric_flag": f"{flag}",
                    "metric_deliverable": True
                    })
            except KeyError:
                job_json["metric"] = [{
                    "metric_name": "purity",
                    "metric_value": f"{value}",
                    "metric_flag": f"{flag}",
                    "metric_deliverable": True
                    }]
            job_json["file"].append({
                "location_uri": f"beluga://{filename}",
                "file_name": f"{os.path.basename(filename)}"
                })
    except FileNotFoundError:
        pass
    if not job_json["file"]:
        job_json = None
    return job_json

def extract_rnaseqc2_rna(sample, prefix_path):
    """
    Extract rnaseqc2 job, files and metrics from the job_output/rnaseqc2/rnaseqc2.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/rnaseqc2/rnaseqc2.*{sample}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/rnaseqc2/rnaseqc2.*{sample}*.o file found")
    job_json = {
        "job_name": f"rnaseqc2.{sample}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": [],
        "metric": []
    }
    try:
        no_expression_profiling_efficiency = no_ribosomal_contamination_count = False
        filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample, "rnaseqc2",  f'{sample}.sorted.mdup.bam.metrics.tsv')
        with open(filename, 'r', encoding="utf-8") as file:
            job_json["job_status"] = "COMPLETED"
            for line in file:
                if line.startswith("Expression Profiling Efficiency"):
                    value = line.split("\t")[1]
                    job_json["metric"].append({
                        "metric_name": "expression_profiling_efficiency",
                        "metric_value": f"{value}",
                        "metric_flag": "PASS",
                        "metric_deliverable": True
                        })
                elif line.startswith("rRNA Rate"):
                    value = line.split("\t")[1]
                    if float(value)>0.35:
                        flag = "FAILED"
                    elif float(value)>0.1:
                        flag = "WARNING"
                    else:
                        flag = "PASS"
                    job_json["metric"].append({
                        "metric_name": "ribosomal_contamination_count",
                        "metric_value": f"{value}",
                        "metric_flag": f"{flag}",
                        "metric_deliverable": True
                        })
            job_json["file"].append({
                "location_uri": f"beluga://{filename}",
                "file_name": f"{os.path.basename(filename)}"
                })
    except FileNotFoundError:
        no_expression_profiling_efficiency = no_ribosomal_contamination_count = True

    if no_expression_profiling_efficiency and no_ribosomal_contamination_count:
        job_json = None

    return job_json

def extract_picard_rna(sample, prefix_path):
    """
    Extract picard job, files and metrics from the job_output/picard_rna_metrics/picard_rna_metrics.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/picard_rna_metrics/picard_rna_metrics.*{sample}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/picard_rna_metrics/picard_rna_metrics.*{sample}*.o file found")
    job_json = {
        "job_name": f"picard.{sample}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": [],
        "metric": []
    }
    try:
        no_bases_over_q30_percent = False
        filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample + '.quality_distribution_metrics')
        tester = re.compile(r'(\d+)\W+(\d+)')
        with open(filename, 'r', encoding="utf-8") as file:
            job_json["job_status"] = "COMPLETED"
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
            value = round((above_30/(above_30+below_30))*100, 2)
            if int(value)<75:
                flag = "FAILED"
            elif int(value)<80:
                flag = "WARNING"
            else:
                flag = "PASS"
            job_json["metric"].append({
                "metric_name": "bases_over_q30_percent",
                "metric_value": f"{value}",
                "metric_flag": f"{flag}",
                "metric_deliverable": True
                })
            job_json["file"].append({
                "location_uri": f"beluga://{filename}",
                "file_name": f"{os.path.basename(filename)}"
                })
    except FileNotFoundError:
        no_bases_over_q30_percent = True

    if no_bases_over_q30_percent:
        job_json = None

    return job_json


def extract_kallisto_rnaseqlight(sample, prefix_path):
    """
    Extract kallisto job, files and metrics from the job_output/kallisto/kallisto.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/kallisto/kallisto.*{sample}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/kallisto/kallisto.*{sample}*.o file found")
    job_json = {
        "job_name": f"kallisto.{sample}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": [],
        "metric": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/kallisto', sample, 'abundance_transcripts.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/kallisto', sample, 'abundance_genes.tsv')
    job_json = add_output_file(filename, job_json, True)
    try:
        filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/kallisto', sample, 'abundance.h5')
        job_json["file"].append({
            "location_uri": f"beluga://{filename}",
            "file_name": f"{os.path.basename(filename)}"
            })
        with h5py.File(filename, 'r') as h5file:
            frequencies = np.asarray(h5file['aux']['fld'])
            values = np.arange(0, len(frequencies), 1)
            # Calculus coming from https://stackoverflow.com/questions/46086663/how-to-get-mean-and-standard-deviation-from-a-frequency-distribution-table
            ord = np.argsort(values)
            cdf = np.cumsum(frequencies[ord])
            median_insert_size = values[ord][np.searchsorted(cdf, cdf[-1] // 2)].astype('float64')
            mean_insert_size = np.around(np.average(values, weights=frequencies), decimals=1)
            job_json["job_status"] = "COMPLETED"
            job_json["metric"].append({
                "metric_name": "mean_insert_size",
                "metric_value": f"{mean_insert_size}",
                "metric_flag": "PASS",
                "metric_deliverable": True
                })
            job_json["metric"].append({
                "metric_name": "median_insert_size",
                "metric_value": f"{median_insert_size}",
                "metric_flag": "PASS",
                "metric_deliverable": True
                })

    except FileNotFoundError:
        logger.warning(f"{filename} not found, no kallisto metrics will be stored.")

    if not job_json["file"]:
        job_json = None

    return job_json

def run_annofuse_rnaseqlight(patient, sample, prefix_path):
    """
    Extract annofuse job and files from the job_output/run_annofuse/run_annoFuse.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/run_annofuse/run_annoFuse.*{sample}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/run_annofuse/run_annoFuse.*{sample}*.o file found")
    job_json = {
        "job_name": f"run_annofuse.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/fusion', sample, 'annoFuse', f'{sample}.putative_driver_fusions.tsv')
    job_json = add_output_file(filename, job_json)

    if not job_json["file"]:
        job_json = None

    return job_json


def tumour_pair_multiqc(patient, prefix_path):
    """
    Extract multiqc job, files and metrics from the job_output/run_pair_multiqc/multiqc.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/run_pair_multiqc/multiqc.*{patient}*.o")), key=os.path.getmtime)[-1]
        with open(latest, 'r', encoding="utf-8") as file:
            job_status = "COMPLETED"
            for line in file:
                if "AccrueTime" in line:
                    job_start = re.findall(r"\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d", line)[0].replace("T", " ")
                elif "EndTime" in line:
                    job_stop = re.findall(r"\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d", line)[-1].replace("T", " ")
    except IndexError:
        logger.warning(f"No job_output/run_pair_multiqc/multiqc.*{patient}*.o file found")
    job_json = {
        "job_name": f"multiqc.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', f'{patient}.multiqc.html')
    if os.path.exists(filename):
        job_json["job_status"] = "COMPLETED"
        job_json["file"].append({
            "location_uri": f"beluga://{filename}",
            "file_name": f"{os.path.basename(filename)}",
            "file_deliverable": True
            })
    else:
        logger.warning(f"No {filename} file found")
        job_json = None

    return job_json

def extract_qualimap(patient, sample, prefix_path):
    """
    Extract qualimap job, files and metrics from the job_output/metrics_dna_sample_qualimap/dna_sample_qualimap.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/metrics_dna_sample_qualimap/dna_sample_qualimap.*{sample}*.o")), key=os.path.getmtime)[-1]
        with open(latest, 'r', encoding="utf-8") as file:
            job_status = "COMPLETED"
            for line in file:
                if "AccrueTime" in line:
                    job_start = re.findall(r"\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d", line)[0].replace("T", " ")
                elif "EndTime" in line:
                    job_stop = re.findall(r"\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d", line)[-1].replace("T", " ")
    except IndexError:
        logger.warning(f"No job_output/metrics_dna_sample_qualimap/dna_sample_qualimap.*{sample}*.o file found")
    job_json = {
        "job_name": f"qualimap.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": [],
        "metric": []
    }

    try:
        no_median_insert_size = no_dedup_mean_coverage = no_aligned_reads_count = False
        filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'qualimap', sample, 'genome_results.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            job_json["file"].append({
                "location_uri": f"beluga://{filename}",
                "file_name": f"{os.path.basename(filename)}"
                })
            job_json["job_status"] = "COMPLETED"
            for line in file:
                # median_insert_size
                if "median insert size" in line:
                    metrics = line.split(" ")
                    value = int(metrics[-1])
                    if value<300:
                        flag = "WARNING"
                    elif value<150:
                        flag = "FAILED"
                    else:
                        flag = "PASS"
                    job_json["metric"].append({
                        "metric_name": "median_insert_size",
                        "metric_value": f"{value}",
                        "metric_flag": f"{flag}",
                        "metric_deliverable": True
                        })
                # dedup_mean_coverage
                elif "mean coverageData" in line:
                    metrics = line.split(" ")
                    value = float(metrics[-1].replace('X', ''))
                    job_json["metric"].append({
                        "metric_name": "dedup_mean_coverage",
                        "metric_value": f"{value}",
                        "metric_flag": "PASS",
                        "metric_deliverable": True
                        })
                elif "number of mapped reads" in line:
                    metrics = line.split(" ")
                    value = int(metrics[-2].replace(",", ""))
                    if value<260000000 and sample.endswith('DN'):
                        flag = "FAILED"
                    elif value<660000000 and sample.endswith('DN'):
                        flag = "WARNING"
                    elif value<530000000 and sample.endswith('DT'):
                        flag = "FAILED"
                    elif value<1330000000 and sample.endswith('DT'):
                        flag = "WARNING"
                    else:
                        flag = "PASS"
                    job_json["metric"].append({
                        "metric_name": "aligned_reads_count",
                        "metric_value": f"{value}",
                        "metric_flag": f"{flag}",
                        "metric_deliverable": True
                        })
    except FileNotFoundError:
        no_median_insert_size = no_dedup_mean_coverage = no_aligned_reads_count = True

    # mean_insert_size
    try:
        no_mean_insert_size = False
        filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'qualimap', sample, "raw_data_qualimapReport", "insert_size_histogram.txt")
        job_json["file"].append({
            "location_uri": f"beluga://{filename}",
            "file_name": f"{os.path.basename(filename)}"
            })
        dataset = np.genfromtxt(fname=filename, delimiter="\t", skip_header=1)
        frequencies = dataset[:, 1]
        values = dataset[:, 0]
        # Calculus coming from https://stackoverflow.com/questions/46086663/how-to-get-mean-and-standard-deviation-from-a-frequency-distribution-table
        value = np.around(np.average(values, weights=frequencies), decimals=1)
        job_json["job_status"] = "COMPLETED"
        job_json["metric"].append({
            "metric_name": "mean_insert_size",
            "metric_value": f"{value}",
            "metric_flag": "PASS",
            "metric_deliverable": True
            })
    except FileNotFoundError:
        no_mean_insert_size = True

    if no_mean_insert_size and no_median_insert_size and no_dedup_mean_coverage and no_aligned_reads_count:
        job_json = None

    return job_json

def extract_picard_tumourpair(patient, sample, prefix_path):
    """
    Extract picard job, files and metrics from the job_output/metrics_dna_picard_metrics/picard_collect_multiple_metrics.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/metrics_dna_picard_metrics/picard_collect_multiple_metrics.*{sample}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/metrics_dna_picard_metrics/picard_collect_multiple_metrics.*{sample}*.o file found")
    try:
        filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'picard_metrics', sample + '.all.metrics.quality_distribution_metrics')
        tester = re.compile(r'(\d+)\W+(\d+)')
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
            value = round((above_30/(above_30+below_30))*100, 2)
            if int(value)<75:
                flag = "FAILED"
            elif int(value)<80:
                flag = "WARNING"
            else:
                flag = "PASS"
            metric_json = {
                "metric_name": "bases_over_q30_percent",
                "metric_value": f"{value}",
                "metric_flag": f"{flag}",
                "metric_deliverable": True
                }
            file_json = {
                "location_uri": f"beluga://{filename}",
                "file_name": f"{os.path.basename(filename)}"
                }
            job_json = {
                "job_name": f"picard.{patient}",
                "job_start": job_start,
                "job_stop": job_stop,
                "job_status": job_status,
                "file": [file_json],
                "metric": [metric_json]
            }
    except FileNotFoundError:
        job_json = None
    return job_json


def gatk_variant_annotator_germline_tumourpair(patient, prefix_path):
    """
    Extract gatk_variant_annotator_germline job and files from the job_output/gatk_variant_annotator_germline/gatk_variant_annotator_germline.others.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/gatk_variant_annotator_germline/gatk_variant_annotator_germline.others.*{patient}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/gatk_variant_annotator_germline/gatk_variant_annotator_germline.others.*{patient}*.o file found")
    job_json = {
        "job_name": f"gatk_variant_annotator_germline.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble', patient, f'{patient}.ensemble.germline.vt.annot.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    if not job_json["file"]:
        job_json = None

    return job_json

def gatk_variant_annotator_somatic_tumourpair(patient, prefix_path):
    """
    Extract gatk_variant_annotator_somatic job and files from the job_output/gatk_variant_annotator_germline/gatk_variant_annotator_somatic.others.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/gatk_variant_annotator_somatic/gatk_variant_annotator_somatic.others.*{patient}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/gatk_variant_annotator_germline/gatk_variant_annotator_somatic.others.*{patient}*.o file found")
    job_json = {
        "job_name": f"gatk_variant_annotator_somatic.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble', patient, f'{patient}.ensemble.somatic.vt.annot.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    if not job_json["file"]:
        job_json = None

    return job_json

def paired_mutect2_tumourpair(patient, prefix_path):
    """
    Extract paired_mutect2 job and files from the job_output/paired_mutect2/gatk_mutect2.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/paired_mutect2/gatk_mutect2.*{patient}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/paired_mutect2/gatk_mutect2.*{patient}*.o file found")
    job_json = {
        "job_name": f"paired_mutect2.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, f'{patient}.mutect2.somatic.vt.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, f'{patient}.mutect2.vt.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    if not job_json["file"]:
        job_json = None

    return job_json

def strelka2_paired_germline_tumourpair(patient, prefix_path):
    """
    Extract strelka2_paired_germline job and files from the job_output/strelka2_paired_germline/strelka2_paired_germline.filter.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/strelka2_paired_germline/strelka2_paired_germline.filter.*{patient}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/strelka2_paired_germline/strelka2_paired_germline.filter.*{patient}*.o file found")
    job_json = {
        "job_name": f"strelka2_paired_germline.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, f'{patient}.strelka2.germline.vt.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, f'{patient}.strelka2.germline.vt.vcf.gz.tbi')
    job_json = add_output_file(filename, job_json)
    if not job_json["file"]:
        job_json = None

    return job_json

def vardict_paired_tumourpair(patient, prefix_path):
    """
    Extract vardict_paired job and files from the job_output/vardict_paired/vardict_paired.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/vardict_paired/vardict_paired.*{patient}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/vardict_paired/vardict_paired.*{patient}*.o file found")
    job_json = {
        "job_name": f"vardict_paired.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, f'{patient}.vardict.germline.vt.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, f'{patient}.vardict.somatic.vt.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    if not job_json["file"]:
        job_json = None

    return job_json

def paired_varscan2_tumourpair(patient, prefix_path):
    """
    Extract paired_varscan2 job and files from the job_output/paired_varscan2/varscan2_somatic.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/paired_varscan2/varscan2_somatic.*{patient}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/paired_varscan2/varscan2_somatic.*{patient}*.o file found")
    job_json = {
        "job_name": f"paired_varscan2.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, f'{patient}.varscan2.germline.vt.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, f'{patient}.varscan2.somatic.vt.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    if not job_json["file"]:
        job_json = None

    return job_json

def cnvkit_batch_tumourpair(patient, prefix_path):
    """
    Extract cnvkit_batch job and files from the job_output/paired_mutect2/gatk_mutect2.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/cnvkit_batch/cnvkit_batch.call.*{patient}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/paired_mutect2/gatk_mutect2.*{patient}*.o file found")
    job_json = {
        "job_name": f"cnvkit_batch.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants', patient, f'{patient}.cnvkit.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    if not job_json["file"]:
        job_json = None

    return job_json

def recalibration_tumourpair(patient, sample, prefix_path):
    """
    Extract recalibration job and files from the job_output/recalibration/gatk_print_reads.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/recalibration/gatk_print_reads.*{sample}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/recalibration/gatk_print_reads.*{sample}*.o file found")
    job_json = {
        "job_name": f"recalibration.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment', sample, f'{sample}.sorted.dup.recal.bam')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment', sample, f'{sample}.sorted.dup.recal.bam.bai')
    job_json = add_output_file(filename, job_json, True)
    if not job_json["file"]:
        job_json = None

    return job_json

def sym_link_final_bam_tumourpair(patient, sample, tumor, prefix_path):
    """
    Extract sym_link_final_bam job and files from the job_output/sym_link_final_bam/sym_link_final_bam.pairs.*{patient}.{t_name}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        if tumor:
            t_name = "Tumor"
        else:
            t_name = "Normal"
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/sym_link_final_bam/sym_link_final_bam.pairs.*{patient}.{t_name}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/sym_link_final_bam/sym_link_final_bam.pairs.*{patient}.{t_name}*.o file found")
    job_json = {
        "job_name": f"sym_link_final_bam.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment', sample, f'{sample}.sorted.dup.recal.bam.md5')
    job_json = add_output_file(filename, job_json)
    if not job_json["file"]:
        job_json = None

    return job_json

def run_pair_multiqc_tumourpair(patient, prefix_path):
    """
    Extract run_pair_multiqc job and files from the job_output/run_pair_multiqc/multiqc.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/run_pair_multiqc/multiqc.*{patient}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/run_pair_multiqc/multiqc.*{patient}*.o file found")
    job_json = {
        "job_name": f"run_pair_multiqc.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', f'{patient}.multiqc.html')
    job_json = add_output_file(filename, job_json, True)
    if not job_json["file"]:
        job_json = None

    return job_json

def report_pcgr_tumourpair(patient, prefix_path):
    """
    Extract report_pcgr job and files from the job_output/report_pcgr/report_pcgr.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/report_pcgr/report_pcgr.*{patient}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/report_pcgr/report_pcgr.*{patient}*.o file found")
    job_json = {
        "job_name": f"report_pcgr.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble', patient, 'pcgr', f'{patient}.pcgr_acmg.grch38.flexdb.html')
    job_json = add_output_file(filename, job_json, True)

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble', patient, 'pcgr', f'{patient}.pcgr_acmg.grch38.maf')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble', patient, 'pcgr', f'{patient}.pcgr_acmg.grch38.snvs_indels.tiers.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble', patient, 'pcgr', f'{patient}.pcgr_acmg.grch38.cna_segments.tsv.gz')
    job_json = add_output_file(filename, job_json, True)
    if not job_json["file"]:
        job_json = None

    return job_json

def parse_o_file(latest):
    """
    Parse the latest job_output file to extract job_status, job_start and job_stop.
    Return: str, str, str
    """
    job_status = None
    job_start = None
    job_stop = None
    with open(latest, 'rb') as file:
        job_status = "COMPLETED"
        file.seek(-4000, 2)
        for line in file.read().decode('utf-8').split("\n"):
            if "AccrueTime" in line:
                job_start = re.findall(r"\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d", line)[0].replace("T", " ")
            elif "EndTime" in line:
                job_stop = re.findall(r"\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d", line)[-1].replace("T", " ")
    return job_status, job_start, job_stop

def add_output_file(filename, job_json, deliverable=False):
    """
    Add output file to the job_json dict if the file exists.
    Return: dict
    """
    if os.path.exists(filename):
        job_json["job_status"] = "COMPLETED"
        job_json["file"].append({
            "location_uri": f"beluga://{filename}",
            "file_name": f"{os.path.basename(filename)}",
            "file_deliverable": deliverable
            })
    return job_json

def dna_bases_over_q30_percent_check(value):
    """ Bases over Q30 (%) DNA metric check """
    if int(value)<75:
        ret = "FAILED"
    elif int(value)<80:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def dna_aligned_reads_count_check(value, tumour):
    """ Aligned Reads Count DNA metric check """
    if int(value)<260000000 and not tumour:
        ret = "FAILED"
    elif int(value)<660000000 and not tumour:
        ret = "WARNING"
    elif int(value)<530000000 and tumour:
        ret = "FAILED"
    elif int(value)<1330000000 and tumour:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def dna_raw_mean_coverage_check(value, tumour):
    """ Mean Coverage DNA metric check """
    if not value:
        ret = "MISSING"
    elif float(value)<30 and not tumour:
        ret = "FAILED"
    elif float(value)<80 and tumour:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def rna_raw_reads_count_check(value):
    """ Clusters RNA metric check """
    if not value:
        ret = "MISSING"
    elif int(value)<80000000:
        ret = "FAILED"
    elif int(value)<100000000:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def dna_raw_duplication_rate_check(value):
    """ Dup. Rate (%) DNA metric check """
    if not value:
        ret = "MISSING"
    elif float(value)>50:
        ret = "FAILED"
    elif float(value)>20:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def median_insert_size_check(value):
    """ Mapped Insert Size (median) metric check """
    if not value:
        ret = "MISSING"
    elif float(value)<300:
        ret = "WARNING"
    elif float(value)<150:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def dna_contamination_check(value):
    """ Contamination (%) DNA metric check """
    if float(value)>5:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def dna_concordance_check(value):
    """ Concordance (%) DNA metric check """
    if float(value)<99:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def dna_tumour_purity_check(value):
    """ Tumour Purity (%) DNA metric check """
    if float(value)<30:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def rna_exonic_rate_check(value):
    """ Exonic Rate (%) RNA metric check """
    if float(value)<0.6:
        ret = "FAILED"
    elif float(value)<0.8:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def rna_ribosomal_contamination_count_check(value):
    """ Ribosomal Contamination (%) RNA metric check """
    if float(value)>0.35:
        ret = "FAILED"
    elif float(value)>0.1:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def rna_ribosomal_contamination_count_compute(rrna_count, rna_aligned_reads_count):
    """ Ribosomal Contamination (%) RNA metric computation """
    return int(rrna_count)/int(rna_aligned_reads_count)

def filter_gatk_rnaseqc2(sample, prefix_path):
    """
    Extract filter_gatk job and files from the job_output/filter_gatk/filter_gatk.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/filter_gatk/filter_gatk.*{sample}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/filter_gatk/filter_gatk.*{sample}*.o file found")
    job_json = {
        "job_name": f"filter_gatk_.{sample}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment',sample,f'{sample}.hc.vt.annot.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment',sample,f'{sample}.hc.vt.annot.flt.vcf.gz')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment',sample,f'{sample}.hc.vt.annot.flt.vcf.gz.tbi')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment',sample,f'{sample}.hc.vt.annot.flt.vcf.gz.md5')
    job_json = add_output_file(filename, job_json, True)

    if not job_json["file"]:
        job_json = None

    return job_json

def recalibration_rnaseqc2(sample, prefix_path):
    """
    Extract recalibration job and files from the job_output/recalibration/gatk_print_reads.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/recalibration/gatk_print_reads.*{sample}*.o")), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/recalibration/gatk_print_reads.*{sample}*.o file found")
    job_json = {
        "job_name": f"gatk_print_reads.{sample}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment/',sample ,f'{sample}.sorted.mdup.split.recal.bam')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment/',sample,f'{sample}.sorted.mdup.split.recal.bam.bai')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment/',sample,f'{sample}.sorted.mdup.split.recal.bam.md5')
    job_json = add_output_file(filename, job_json, True)

    if not job_json["file"]:
        job_json = None

    return job_json

def report_pcgr_rnaseqc2(sample, prefix_path):
    """
    Extract report_pcgr job and files from the job_output/report_pcgr/report_pcgr.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/report_pcgr/report_pcgr.*{sample}*.o" )), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/report_pcgr/report_pcgr.*{sample}*.o file found")
    job_json = {
        "job_name": f"report_pcgr.{sample}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment/', sample, 'pcgr', f'{sample}.pcgr_acmg.grch38.flexdb.html')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment/', sample, 'pcgr', f'{sample}.pcgr_acmg.grch38.maf')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/alignment/', sample, 'pcgr', f'{sample}.pcgr_acmg.grch38.snvs_indels.tiers.tsv')
    job_json = add_output_file(filename, job_json, True)

    if not job_json["file"]:
        job_json = None

    return job_json


def rnaseqc2_multiqc(sample, prefix_path):
    """
    Extract run_pair_multiqc job and files from the job_output/run_pair_multiqc/multiqc.*{sample}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/run_pair_multiqc/multiqc.*{sample}*.o")), key=os.path.getmtime)[-1]
        with open(latest, 'r', encoding="utf-8") as file:
            job_status = "COMPLETED"
            for line in file:
                if "AccrueTime" in line:
                    job_start = re.findall(r"\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d", line)[0].replace("T", " ")
                elif "EndTime" in line:
                    job_stop = re.findall(r"\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d", line)[-1].replace("T", " ")
    except IndexError:
        logger.warning(f"No job_output/multiqc/multiqc.*{sample}*.o file found")
    job_json = {
        "job_name": f"multiqc.{sample}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/metrics/multiqc_by_sample', sample, f'{sample}.multiqc.html')
    if os.path.exists(filename):
        job_json["job_status"] = "COMPLETED"
        job_json["file"].append({
            "location_uri": f"beluga://{filename}",
            "file_name": f"{os.path.basename(filename)}",
            "file_deliverable": True
            })
    else:
        logger.warning(f"No {filename} file found")
        job_json = None

    return job_json

def gridss_paired_somatic(patient, sample, prefix_path):
    """
    Extract gridss_paired_somatic job and files from the job_output/gridss_paired_somatic/gridss_paired_somatic.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/gridss_paired_somatic/gridss_paired_somatic.*{patient}*.o" )), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/gridss_paired_somatic/gridss_paired_somatic.*{patient}*.o file found")
    job_json = {
        "job_name": f"gridss_paired.somatic.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'gridss', f'{sample}.gridss.vcf.gz')
    job_json = add_output_file(filename, job_json, True)

    if not job_json["file"]:
        job_json = None

    return job_json

def gripss_filter_somatic(patient, sample, prefix_path):
    """
    Extract gripss_filter_somatic job and files from the job_output/gridss_paired_somatic/gripss_filter.somatic.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/gridss_paired_somatic/gripss_filter.somatic.*{patient}*.o" )), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/gridss_paired_somatic/gripss_filter.somatic.*{patient}*.o file found")
    job_json = {
        "job_name": f"gripss_filter.somatic.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'gridss', f'{sample}.gripss.filtered.somatic.vcf.gz')
    job_json = add_output_file(filename, job_json, True)

    if not job_json["file"]:
        job_json = None

    return job_json

def purple_purity_sv(patient, sample, prefix_path):
    """
    Extract purple_purity_sv job and files from the job_output/purple_sv/purple.purity.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/purple_sv/purple.purity.*{patient}*.o" )), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/purple_sv/purple.purity.*{patient}*.o file found")
    job_json = {
        "job_name": f"purple.purity.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'purple/plot', f'{sample}.circos.png')
    job_json = add_output_file(filename, job_json, True)

    if not job_json["file"]:
        job_json = None

    return job_json

def gripss_filter_germline(patient,sample,prefix_path):
    """
    Extract gripss_filter_germline job and files from the job_output/gridss_paired_somatic/gripss_filter.germline.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/gridss_paired_somatic/gripss_filter.germline.*{patient}*.o" )), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/gridss_paired_somatic/gripss_filter.germline.*{patient}*.o file found")
    job_json = {
        "job_name": f"gripss_filter.germline.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'gridss', f'{sample}.gripss.filtered.germline.vcf.gz')
    job_json = add_output_file(filename, job_json, True)

    if not job_json["file"]:
        job_json = None

    return job_json

def linx_annotations_somatic(patient, sample, prefix_path):
    """
    Extract linx_annotations_somatic job and files from the job_output/linx_annotations_somatic/linx_annotations_somatic.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/linx_annotations_somatic/linx_annotations_somatic.*{patient}*.o" )), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/linx_annotations_somatic/linx_annotations_somatic.*{patient}*.o file found")
    job_json = {
        "job_name": f"linx_annotations_somatic.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.breakend.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.clusters.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.driver.catalog.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.drivers.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.fusion.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.links.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.svs.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.vis_copy_number.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.vis_fusion.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.vis_gene_exon.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.vis_protein_domain.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.vis_segments.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.vis_sv_data.tsv')
    job_json = add_output_file(filename, job_json, True)

    if not job_json["file"]:
        job_json = None

    return job_json

def linx_annotations_germline(patient, sample, prefix_path):
    """
    Extract linx_annotations_germline job and files from the job_output/linx_annotations_germline/linx_annotations_germline.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/linx_annotations_germline/linx_annotations_germline.*{patient}*.o" )), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/linx_annotations_germline/linx_annotations_germline.*{patient}*.o file found")
    job_json = {
        "job_name": f"linx_annotations_germline.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.germline.clusters.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.germline.disruption.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.germline.driver.catalog.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.germline.links.tsv')
    job_json = add_output_file(filename, job_json, True)
    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx', f'{sample}.linx.germline.svs.tsv')
    job_json = add_output_file(filename, job_json, True)


    if not job_json["file"]:
        job_json = None

    return job_json

def linx_plot(patient, prefix_path):
    """
    Extract linx_plot job and files from the job_output/linx_plot/linx_plot.*{patient}*.o file.
    Return: dict
    """
    job_status = None
    job_start = None
    job_stop = None
    try:
        latest = sorted(glob.glob(os.path.join(prefix_path, f"C3G/projects/MOH_PROCESSING/MAIN/job_output/linx_plot/linx_plot.*{patient}*.o" )), key=os.path.getmtime)[-1]
        job_status, job_start, job_stop = parse_o_file(latest)
    except IndexError:
        logger.warning(f"No job_output/linx_plot/plinx_plot.*{patient}*.o file found")
    job_json = {
        "job_name": f"linx.plot.{patient}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": job_status,
        "file": []
    }

    filename = os.path.join(prefix_path, 'C3G/projects/MOH_PROCESSING/MAIN/SVariants/', patient, 'linx/plot')
    job_json = add_output_file(filename, job_json, True)

    if not job_json["file"]:
        job_json = None

    return job_json

if __name__ == '__main__':
    main()
