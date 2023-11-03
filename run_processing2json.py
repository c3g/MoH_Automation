#!/usr/bin/env python3

import argparse
import json
import glob
import sys
import re
import csv
import os
import hashlib
import logging
from datetime import datetime

logging.basicConfig(format='%(levelname)s: %(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(prog='run_processing2json.py', description="Creates json file for project tracking database for a given run processing.")
    parser.add_argument('--input', required=True, help="Input align_bwa_mem.csv file from Run Processing.")
    parser.add_argument('--output', required=False, help="Output json filename (Default: <input_filename>.json).")
    args = parser.parse_args()

    if not args.output:
        output = f"{os.path.basename(args.input).split('.')[0]}.json"
    else:
        output = args.output

    jsonify_run_processing(args.input, output)

def jsonify_run_processing(input_csv, output):
    readset_dict = {}
    sample_dict = {}
    run_list = []
    with open(input_csv, 'rt') as run_file_in:
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
            cohort = result.group(2)
            institution = result.group(3)
            patient_json = {
                "patient_fms_id": None,
                "patient_name": patient,
                "patient_cohort": cohort,
                "patient_institution": institution,
                "sample": []
                }
            if sample.endswith("T"):
                sample_tumour = True
            else:
                sample_tumour = False
            sample_json = {
                "sample_fms_id": None,
                "sample_name": sample,
                "sample_tumour": sample_tumour,
                "readset": []
                }

            copylist = os.path.join(os.path.dirname(input_csv), f"{os.path.basename(input_csv).split('.')[0]}.copylist.txt")
            if not os.path.isfile(copylist):
                raise Exception(f"File {copylist} not found; required to fing raw data (bams/fastqs) location")
            fastq1 = fastq2 = ""
            with open(copylist, 'r') as file:
                for line in file:
                    if sample in line:
                        fields = line.split(",")
                        file_path = fields[3].strip()
                        if file_path.endswith(".bam"):
                            file_json = [
                                {
                                    "location_uri": f"abacus://{file_path}",
                                    "file_name": f"{os.path.basename(file_path)}",
                                    "file_deliverable": True
                                    }
                                ]
                            break
                        elif file_path.endswith(".fastq.gz"):
                            if "R1" in file_path:
                                fastq1 = os.path.basename(file_path)
                                fastq1_location_uri = file_path
                            elif "R2" in file_path:
                                fastq2 = os.path.basename(file_path)
                                fastq2_location_uri = file_path
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
            patient_json["sample"].append(sample_json)
            json_output["patient"].append(patient_json)

    with open(output, 'w', encoding='utf-8') as file:
        json.dump(json_output, file, ensure_ascii=False, indent=4)

    return readset_dict, sample_dict


def dna_raw_mean_coverage_check(value, tumour):
    if not value:
        raise Exception(f"Missing 'Mean Coverage': {value}")
    if float(value)<30 and not tumour:
        ret = "FAILED"
    elif float(value)<80 and tumour:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def rna_raw_reads_count_check(value):
    if not value:
        raise Exception(f"Missing 'Clusters': {value}")
    if int(value)<80000000:
        ret = "FAILED"
    elif int(value)<100000000:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def dna_raw_duplication_rate_check(value):
    if not value:
        raise Exception(f"Missing 'Dup. Rate (%)': {value}")
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
    if not value:
        raise Exception(f"Missing 'Mapped Insert Size (median)': {value}")
    if float(value)<300:
        ret = "WARNING"
    elif float(value)<150:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

if __name__ == '__main__':
    main()
