#!/usr/bin/env python3

import argparse
import json
import re
import csv
import os
import logging
from datetime import datetime

logging.basicConfig(format='%(levelname)s: %(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    """ Main """
    parser = argparse.ArgumentParser(prog='run_processing2json.py', description="Creates json file for project tracking database for a given run processing.")
    parser.add_argument('-i', '--input', required=True, help="Input align_bwa_mem.csv file from Run Processing.")
    parser.add_argument('-o', '--output', required=False, help="Output json filename (Default: <input_filename>.json).")
    parser.add_argument('-l', '--lane', required=False, help="Only considers lane(s) provided for json creation.", nargs='+')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--sample', required=False, help="Only considers sample(s) provided for json creation.", nargs='+')
    group.add_argument('-x', '--xsample', required=False, help="Ignores sample(s) provided for json creation.", nargs='+')
    args = parser.parse_args()

    if not args.output:
        output = f"{os.path.basename(args.input).split('.')[0]}.json"
    else:
        output = args.output
    if args.lane:
        lanes = list(args.lane)
    else:
        lanes = ["1", "2", "3", "4"]

    samples = []
    input_csv = args.input
    run_list = []
    with open(input_csv, 'rt') as run_file_in:
        reader = csv.DictReader(run_file_in)
        for row in reader:
            run_list.append(row)
            samples.append(row['Sample Name'])

    if args.sample:
        samples = list(args.sample)
    elif args.xsample:
        samples = list(set(samples).difference(list(args.xsample)))

    jsonify_run_processing(input_csv, run_list, output, lanes, samples)

def jsonify_run_processing(input_csv, run_list, output, lanes, samples):
    """ Writing RUn Processing json based on csv"""
    readset_dict = {}
    sample_dict = {}
    # change fms_id to ext_id and ext_src once Database is updated in prod
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
        if sample.startswith("MoHQ") and run_row['Lane'] in lanes and sample in samples:
            result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX|IQ|HM)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
            patient = result.group(1)
            cohort = result.group(2)
            institution = result.group(3)
            # change fms_id to ext_id and ext_src once Database is updated in prod
            patient_json = {
                "patient_fms_id": None,
                "patient_name": patient,
                "patient_cohort": cohort,
                "patient_institution": institution,
                "sample": []
                }
            # if sample.endswith("T"):
            #     sample_tumour = True
            # else:
            #     sample_tumour = False
            sample_tumour = sample.endswith("T")
            # change fms_id to ext_id and ext_src once Database is updated in prod
            sample_json = {
                "sample_fms_id": None,
                "sample_name": sample,
                "sample_tumour": sample_tumour,
                "readset": []
                }

            copylist = os.path.join(os.path.dirname(input_csv), f"{os.path.basename(input_csv).split('.')[0]}.copylist.txt")
            if not os.path.isfile(copylist):
                raise Exception(f"File {copylist} not found; required to find raw data (bams/bais/fastqs) location")
            fastq1 = fastq2 = bam = bai = ""
            with open(copylist, 'r') as file:
                for line in file:
                    if sample in line:
                        fields = line.split(",")
                        file_path = fields[3].strip()
                        if re.search(r'ba.$', file_path):
                            if file_path.endswith(".bam"):
                                bam = os.path.basename(file_path)
                                bam_location_uri = file_path
                            elif file_path.endswith(".bai"):
                                bai = os.path.basename(file_path)
                                bai_location_uri = file_path
                            if bam and bai:
                                file_json = [
                                    {
                                        "location_uri": f"abacus://{bam_location_uri}",
                                        "file_name": f"{bam}",
                                        "file_deliverable": True
                                        },
                                    {
                                        "location_uri": f"abacus://{bai_location_uri}",
                                        "file_name": f"{bai}",
                                        "file_deliverable": True
                                        }
                                    ]
                                break
                        elif file_path.endswith(".fastq.gz"):
                            if "_R1_" in file_path:
                                fastq1 = os.path.basename(file_path)
                                fastq1_location_uri = file_path
                            elif "_R2_" in file_path:
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
                raw_reads_count_flag = rna_raw_reads_count_check(sample, run_row['Clusters'])
            raw_duplication_rate_flag = "PASS"
            if run_row['Library Type'] != "RNASeq":
                raw_duplication_rate_flag = dna_raw_duplication_rate_check(sample, run_row['Dup. Rate (%)'])
            raw_median_insert_size_flag = median_insert_size_check(sample, run_row['Mapped Insert Size (median)'])
            raw_mean_insert_size_flag = "PASS"
            raw_mean_coverage_flag = "PASS"
            if run_row['Library Type'] != "RNASeq":
                raw_mean_coverage_flag = dna_raw_mean_coverage_check(sample, run_row['Mean Coverage'], sample_tumour)
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
                "experiment_nucleic_acid_type": "RNA" if run_row['Library Type'] == "RNASeq" else "DNA",
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


def dna_raw_mean_coverage_check(sample, value, tumour):
    """ Mean Coverage DNA metric check """
    if not value:
        raise Exception(f"Missing 'Mean Coverage' value for {sample}")
    if float(value)<30 and not tumour:
        ret = "FAILED"
    elif float(value)<80 and tumour:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def rna_raw_reads_count_check(sample, value):
    """ Clusters RNA metric check """
    if not value:
        raise Exception(f"Missing 'Clusters' value for {sample}")
    if int(value)<80000000:
        ret = "FAILED"
    elif int(value)<100000000:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def dna_raw_duplication_rate_check(sample, value):
    """ Dup. Rate (%) DNA metric check """
    if not value:
        raise Exception(f"Missing 'Dup. Rate (%)' value for {sample}")
    if not value:
        ret = "FAILED"
    elif float(value)>50:
        ret = "FAILED"
    elif float(value)>20:
        ret = "WARNING"
    else:
        ret = "PASS"
    return ret

def median_insert_size_check(sample, value):
    """ Mapped Insert Size (median) metric check """
    if not value:
        raise Exception(f"Missing 'Mapped Insert Size (median)' value for {sample}")
    if float(value)<300:
        ret = "WARNING"
    elif float(value)<150:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

if __name__ == '__main__':
    main()
