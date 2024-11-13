#!/usr/bin/env python3

import argparse
import json
import re
import glob
import os
import logging
from datetime import datetime

logging.basicConfig(format='%(levelname)s: %(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        prog='run_processing2json.py',
        description="Creates json file for project tracking database from a Run generated by GenPipes Run Processing and Freezeman."
        )
    parser.add_argument(
        '-i',
        '--input',
        required=True,
        help="Folder of the Run to be converted. To be formatted like this /lb/robot/research/freezeman-processing/<sequencer>/<year>/<run_folder>."
        )
    parser.add_argument(
        '-o',
        '--output',
        default=None,
        help="Output json filename (Default: .json)."
        )
    parser.add_argument(
        '-l',
        '--lane',
        nargs='+',
        default=["1", "2", "3", "4", "5", "6", "7", "8"],
        help="Only considers lane(s) provided for json creation."
        )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '-s',
        '--sample',
        nargs='+',
        help="Only considers sample(s) provided for json creation."
        )
    group.add_argument(
        '-x',
        '--xsample',
        nargs='+',
        help="Ignores sample(s) provided for json creation."
        )
    return parser.parse_args()

def load_json_file(file_path):
    """Load a JSON file."""
    try:
        with open(file_path, mode="r") as file:
            return json.load(file)
    except FileNotFoundError:
        logger.error(f"JSON file not found: {file_path}")
    except json.JSONDecodeError:
        logger.error(f"Error decoding JSON file: {file_path}")
    return None

def parse_sample_name(sample_name):
    """Parse the sample name to extract specimen, cohort, and institution."""
    result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX|HM)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample_name)
    return result.group(1), result.group(2), result.group(3)

def jsonify_run_processing(input_run_folder, fms_json, lanes_json, output, lanes, samples):
    """
    Converts the Run Processing csv file to a json file for project tracking database.
    Args:
        input_run_folder (str): Path to the Run Processing csv file.
        fms_json (json): JSON file from Freezeman.
        output (str): Path to the output json file.
        lanes (list): List of lanes to consider.
        samples (list): List of samples to consider.
    """
    readset_dict = {}
    sample_dict = {}
    json_output = {
            "operation_platform": "abacus",
            "project_ext_id": None,
            "project_ext_src": None,
            "project_name": "MOH-Q",
            "run_ext_id": fms_json["run_obj_id"],
            "run_ext_src": "FREEZEMAN",
            "run_name": input_run_folder.split("/")[-1],
            "run_instrument": "novaseq",
            "run_date": f"{datetime.strptime(fms_json['run_start_date'], '%Y-%m-%d')}",
            "specimen": []
            }

    for _, lane_json in lanes_json.items():
        for readset_key, readset in lane_json["readsets"].items():
            sample_name = readset["sample_name"]
            if sample_name.startswith("MoHQ") and lane_json['lane'] in lanes and sample_name in samples:
                specimen, cohort, institution = parse_sample_name(sample_name)
                # Check if the specimen is already in json_output["specimen"]
                specimen_names = [spec["specimen_name"] for spec in json_output["specimen"]]
                if specimen in specimen_names:
                    # Specimen is present, find its position
                    position = specimen_names.index(specimen)
                    specimen_json = json_output["specimen"][position]
                else:
                    # Specimen is not present, add it to json_output["specimen"]
                    specimen_json = {
                        "specimen_ext_id": None,
                        "specimen_ext_src": None,
                        "specimen_name": specimen,
                        "specimen_cohort": cohort,
                        "specimen_institution": institution,
                        "sample": []
                        }
                    json_output["specimen"].append(specimen_json)

                # Check if the sample is already in specimen_json["sample"]
                sample_names = [spec["sample_name"] for spec in specimen_json["sample"]]
                if sample_name in sample_names:
                    # sample is present, find its position
                    position = sample_names.index(sample_name)
                    sample_json = specimen_json["sample"][position]
                else:
                    # sample is not present, add it to specimen_json["sample"]
                    sample_json = {
                        "sample_ext_id": int(readset["derived_sample_obj_id"]),
                        "sample_ext_src": "FREEZEMAN",
                        "sample_name": sample_name,
                        "sample_tumour": sample_name.endswith("T"),
                        "readset": []
                        }
                    specimen_json["sample"].append(sample_json)
                if sample_name.endswith("RT"):
                    fastq1 = readset["fastq_1"]["final_path"]
                    fastq2 = readset["fastq_2"]["final_path"]
                    file_json = [
                        {
                            "location_uri": f"abacus://{fastq1}",
                            "file_name": f"{os.path.basename(fastq1)}",
                            "file_extra_metadata": {"read_type": "R1"},
                            "file_deliverable": True
                            },
                        {
                            "location_uri": f"abacus://{fastq1}",
                            "file_name": f"{os.path.basename(fastq1)}",
                            "file_extra_metadata": {"read_type": "R2"},
                            "file_deliverable": True
                            }
                        ]
                else:
                    bam = readset["bam"]["final_path"]
                    bai = readset["bai"]["final_path"]
                    file_json = [
                        {
                            "location_uri": f"abacus://{bam}",
                            "file_name": f"{os.path.basename(bam)}",
                            "file_deliverable": True
                            },
                        {
                            "location_uri": f"abacus://{bai}",
                            "file_name": f"{os.path.basename(bai)}",
                            "file_deliverable": True
                            }
                        ]
                for run_v in lane_json["run_validation"]:
                    if run_v.get("sample") == readset_key:
                        raw_reads_count = run_v["qc"]["nb_reads"]
                        if raw_reads_count == 0:
                            raw_reads_count_flag = "FAILED"
                        else:
                            raw_reads_count_flag = "PASS"
                        if readset["library_type"] == "RNASeq":
                            raw_reads_count_flag = rna_raw_reads_count_check(sample_name, raw_reads_count)
                        try:
                            raw_duplication_rate = run_v["qc"]["duplication_rate"]
                            raw_duplication_rate_flag = "PASS"
                        except KeyError:
                            raw_duplication_rate = None
                            raw_duplication_rate_flag = "MISSING"
                        if readset["library_type"] != "RNASeq":
                            raw_duplication_rate_flag = dna_raw_duplication_rate_check(sample_name, raw_duplication_rate)
                        try:
                            raw_median_insert_size = run_v["alignment"]["median_aligned_insert_size"]
                            raw_median_insert_size_flag = median_insert_size_check(sample_name, raw_median_insert_size)
                        except KeyError:
                            raw_median_insert_size = None
                            raw_median_insert_size_flag = "MISSING"
                        raw_mean_insert_size = run_v["alignment"]["average_aligned_insert_size"]
                        raw_mean_insert_size_flag = "PASS"
                        raw_mean_coverage = run_v["alignment"]["mean_coverage"]
                        if not raw_mean_coverage:
                            raw_mean_coverage_flag = "MISSING"
                        else:
                            raw_mean_coverage_flag = "PASS"
                            if readset["library_type"] != "RNASeq":
                                raw_mean_coverage_flag = dna_raw_mean_coverage_check(sample_name, raw_mean_coverage, sample_name.endswith("T"))
                        metric_json = [
                            {
                                "metric_name": "raw_reads_count",
                                "metric_value": raw_reads_count,
                                "metric_flag": raw_reads_count_flag,
                                "metric_deliverable": True
                                },
                            {
                                "metric_name": "raw_duplication_rate",
                                "metric_value": raw_duplication_rate,
                                "metric_flag": raw_duplication_rate_flag
                                },
                            {
                                "metric_name": "raw_median_insert_size",
                                "metric_value": raw_median_insert_size,
                                "metric_flag": raw_median_insert_size_flag
                                },
                            {
                                "metric_name": "raw_mean_insert_size",
                                "metric_value": raw_mean_insert_size,
                                "metric_flag": raw_mean_insert_size_flag
                                },
                            {
                                "metric_name": "raw_mean_coverage",
                                "metric_value": raw_mean_coverage,
                                "metric_flag": raw_mean_coverage_flag
                                }
                            ]

                readset_name = f"{sample_name}.{lane_json['run']}_{lane_json['lane']}"
                readset_dict[readset_name] = (specimen, sample_name)
                # Check if the readset is already in sample_json["readset"]
                readset_names = [spec["readset_name"] for spec in sample_json["readset"]]
                if readset_name in readset_names:
                    print(f"Duplicate readset: {readset_name}")
                else:
                    # readset is not present, add it to specimen_json["readset"]
                    for fms_s in fms_json["samples"]:
                        if fms_s.get("sample_name") == sample_name:
                            library_kit = fms_s["library_kit"]
                    readset_json = {
                        "experiment_sequencing_technology": None,
                        "experiment_type": f"{readset["library_type"]}",
                        "experiment_nucleic_acid_type": "RNA" if readset["library_type"] == "RNASeq" else "DNA",
                        "experiment_library_kit": library_kit,
                        "experiment_kit_expiration_date": None,
                        "readset_name": readset_name,
                        "readset_lane": f"{lane_json['lane']}",
                        "readset_adapter1": f"{readset['barcodes'][0]['ADAPTERi7']}",
                        "readset_adapter2": f"{readset['barcodes'][0]['ADAPTERi5']}",
                        "readset_sequencing_type": f"{lane_json['sequencing_method']}",
                        "readset_quality_offset": "33",
                        "file": file_json,
                        "metric": metric_json
                        }
                    sample_json["readset"].append(readset_json)

    with open(output, 'w', encoding='utf-8') as file:
        json.dump(json_output, file, ensure_ascii=False, indent=4)

    return readset_dict, sample_dict


def dna_raw_mean_coverage_check(sample, value, tumour):
    """ Mean Coverage DNA metric check """
    if not value:
        ret = "MISSING"
        logger.warning(f"Missing 'mean_coverage' value for {sample} from json.")
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
        ret = "MISSING"
        logger.warning(f"Missing 'nb_reads' value for {sample} from json.")
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
        ret = "MISSING"
        logger.warning(f"Missing 'duplication_rate' value for {sample} from json.")
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
        ret = "MISSING"
        logger.warning(f"Missing 'median_aligned_insert_size' value for {sample} from json.")
    if float(value)<300:
        ret = "WARNING"
    elif float(value)<150:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def main():
    """ Main """
    args = parse_arguments()

    output = args.output or f"{os.path.basename(args.input)}.json"
    lanes = args.lane

    fms_json = load_json_file(glob.glob(os.path.join(args.input, "*.json"))[0])
    if not fms_json:
        return

    lanes_json_files = glob.glob(os.path.join(args.input, "report", "*.run_validation_report.json"))
    lanes_json = {os.path.basename(file): load_json_file(file) for file in lanes_json_files}

    if args.sample:
        samples = list(args.sample)
    elif args.xsample:
        all_samples = [sample["sample_name"] for sample in fms_json["samples"]]
        samples = list(set(all_samples).difference(set(args.xsample)))
    else:
        samples = [sample["sample_name"] for sample in fms_json["samples"]]

    jsonify_run_processing(args.input, fms_json, lanes_json, output, lanes, samples)

if __name__ == '__main__':
    main()
