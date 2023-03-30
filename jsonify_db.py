#!/usr/bin/env python3

import json
import glob
import sys
import re
import csv
import os
from datetime import datetime

def main():
    main_raw_reads_folder = "/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/raw_reads"
    patient_dict = get_patient_dict(main_raw_reads_folder)
    readset_dict, sample_dict = jsonify_run_processing(patient_dict)
    jsonify_transfer(sample_dict)


def jsonify_run_processing(patient_dict):
    readset_dict = {}
    sample_dict = {}
    run_metrics_csv = "/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/*.csv"
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
                result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
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
                        # print(json_output["patient"][existing_patient_pos])
                        # print(existing_patient, "\n")
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
                        if patient_dict[patient][sample][2].endswith("T"):
                            sample_tumour = True
                        else:
                            sample_tumour = False
                        sample_json = {
                            "sample_fms_id": None,
                            "sample_name": sample_name,
                            "sample_tumour": sample_tumour,
                            "readset": []
                            }

                    # transfer_folder = '/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*'
                    transfer_folder = 'transfer/*'
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
                                                "file_name": f"{os.path.basename(fields[0])}"
                                                }
                                            ]
                                        break
                                    elif ".fastq" in line:
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
                                                "file_extra_metadata": {"read_type": "R1"}
                                                },
                                            {
                                                "location_uri": f"abacus://{fastq2_location_uri}",
                                                "file_name": f"{fastq2}",
                                                "file_extra_metadata": {"read_type": "R2"}
                                                }
                                            ]
                                            break
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
                            "metric_flag": raw_reads_count_flag
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
                        sample_dict[sample_name].append(readset_name)
                    else:
                        sample_dict[sample_name] = [readset_name]
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

def jsonify_transfer(sample_dict):
    # transfer_folder = '/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*'
    transfer_folder = 'transfer/*'
    # transfer_dict = {}
    for filename in glob.glob(transfer_folder):
        transfer_dict = {}
        json_output = {
            "operation_platform": "beluga",
            "operation_cmd_line": f"globus transfer --submission-id $sub_id --label $label --batch /lb/project/mugqic/projects/MOH/TEMP/{os.path.basename(filename)} 6c66d53d-a79d-11e8-96fa-0a6d4e044368 278b9bfe-24da-11e9-9fa2-0a06afd4a22e",
            "readset": []
            }
        with open(filename, 'r') as file:
            # print(filename)
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
            if "fastq1_src_location_uri" in transfer_dict[readset]:
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
            readset_json = {
                "readset_name": readset,
                "file": file_json,
                }
            json_output["readset"].append(readset_json)
            # print(json.dumps(readset_json, indent=4))
        if transfer_dict:
            with open(f"jsons/{os.path.basename(filename).replace('.txt', '.json')}", 'w', encoding='utf-8') as f:
                json.dump(json_output, f, ensure_ascii=False, indent=4)


def jsonify_genpipes():
    pass


def get_patient_dict(main_raw_reads_folder):
    patient_dict = {}
    for sample in os.listdir(main_raw_reads_folder):
        if sample.startswith("MoHQ"):
            result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
            patient = result.group(1)
            cohort = result.group(2)
            institution = result.group(3)
            sample_type = sample[-2:]
            try:
                patient_dict[patient][sample] = (patient, sample, sample_type, cohort, institution)
            except KeyError:
                patient_dict[patient] = {sample: (patient, sample, sample_type, cohort, institution)}
    return patient_dict


def dna_bases_over_q30_percent_check(value):
    if int(value)<75:
        ret = "FAILED"
    elif int(value)<80:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def dna_aligned_reads_count_check(value, tumour):
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
    if float(value)<30 and not tumour:
        ret = "FAILED"
    elif float(value)<80 and tumour:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def rna_raw_reads_count_check(value):
    if int(value)<80000000:
        ret = "FAILED"
    elif int(value)<100000000:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def dna_raw_duplication_rate_check(value):
    if not value:
        ret = "FAILED"
    elif float(value)>50:
        ret = "FAILED"
    elif float(value)>20:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def median_insert_size_check(value):
    if float(value)<300:
        ret = "WARNING"
    elif float(value)<150:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def dna_contamination_check(value):
    if float(value)>5:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def dna_concordance_check(value):
    if float(value)<99:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def dna_tumour_purity_check(value):
    if float(value)<30:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def rna_exonic_rate_check(value):
    if float(value)<0.6:
        ret = "FAILED"
    elif float(value)<0.8:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def rna_ribosomal_contamination_count_check(value):
    if float(value)>0.35:
        ret = "FAILED"
    elif float(value)>0.1:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def rna_ribosomal_contamination_count_compute(rrna_count, rna_aligned_reads_count):
    return int(rrna_count)/int(rna_aligned_reads_count)

if __name__ == '__main__':
    main()
