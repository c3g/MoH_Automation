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

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        prog='mgc_mock_run_processing2json.py',
        description="Creates json file for project tracking database from a Mock Run from Genome Quebec."
        )
    parser.add_argument(
        '-i',
        '--input',
        required=True,
        help="Input align_bwa_mem.csv file from Run Processing."
        )
    parser.add_argument(
        '-o',
        '--output',
        default=None,
        help="Output json filename (Default: .json)."
        )
    # parser.add_argument(
    #     '-n',
    #     '--nucleic_acid_type',
    #     choices=['DNA', 'RNA'],
    #     default='ALL',
    #     help="Nucleic Acid type, either DNA or RNA (Default: ALL)."
    #     )
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

def jsonify_run_processing(input_csv, run_list, output, lanes, samples):
    """ Writing RUn Processing json based on csv"""
    readset_dict = {}
    sample_dict = {}
    json_output = {
            "operation_platform": "abacus",
            "project_ext_id": None,
            "project_ext_src": None,
            "project_name": "MOH-Q",
            "run_ext_id": None,
            "run_ext_src": None,
            "run_name": f"{run_list[0]['Processing Folder Name']}",
            "run_instrument": "novaseq",
            "run_date": f"{datetime.strptime(run_list[0]['Processing Folder Name'][0:6], '%y%m%d')}",
            "specimen": []
            }
    for run_row in run_list:
        sample = run_row['Sample Name']
        if sample.startswith("MoHQ") and run_row['Lane'] in lanes and sample in samples:
            result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX|HM|CQ|IQ)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
            specimen = result.group(1)
            cohort = result.group(2)
            institution = result.group(3)
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

            sample_tumour = sample.endswith("T")
            # Check if the sample is already in specimen_json["sample"]
            sample_names = [spec["sample_name"] for spec in specimen_json["sample"]]
            if sample in sample_names:
                # sample is present, find its position
                position = sample_names.index(sample)
                sample_json = specimen_json["sample"][position]
            else:
                # sample is not present, add it to specimen_json["sample"]
                sample_json = {
                    "sample_ext_id": None,
                    "sample_ext_src": None,
                    "sample_name": sample,
                    "sample_tumour": sample_tumour,
                    "readset": []
                    }
                specimen_json["sample"].append(sample_json)

            copylist = os.path.join(os.path.dirname(input_csv), f"{os.path.basename(input_csv).split('.')[0]}.copylist.txt")
            if not os.path.isfile(copylist):
                raise Exception(f"File {copylist} not found; required to find raw data (bams/bais/fastqs) location")
            fastq1 = fastq2 = bam = bai = ""
            with open(copylist, 'r') as file:
                for line in file:
                    if re.search(fr"{sample}/run{run_row['Run ID']}_{run_row['Lane']}.*\.ba(m|i)$", line):
                        fields = line.split(",")
                        file_path = fields[3].strip()
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
                                    "file_md5sum": compute_md5(bam),
                                    "file_deliverable": True
                                    },
                                {
                                    "location_uri": f"abacus://{bai_location_uri}",
                                    "file_name": f"{bai}",
                                    "file_md5sum": compute_md5(bai),
                                    "file_deliverable": True
                                    }
                                ]
                            break
                    elif re.search(fr"Unaligned\.{run_row['Lane']}/.*/Sample_{sample}.*\.fastq\.gz$", line):
                        fields = line.split(",")
                        file_path = fields[3].strip()
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
                                "file_md5sum": compute_md5(fastq1),
                                "file_extra_metadata": {"read_type": "R1"},
                                "file_deliverable": True
                                },
                            {
                                "location_uri": f"abacus://{fastq2_location_uri}",
                                "file_name": f"{fastq2}",
                                "file_md5sum": compute_md5(fastq2),
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
            readset_name = f"{sample}.{run_row['Run ID']}_{run_row['Lane']}"
            readset_dict[readset_name] = (specimen, sample)
            # Check if the readset is already in sample_json["readset"]
            readset_names = [spec["readset_name"] for spec in sample_json["readset"]]
            if readset_name in readset_names:
                print(f"Duplicate readset: {readset_name}")
            else:
                # readset is not present, add it to specimen_json["readset"]
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

            # sample_json["readset"].append(readset_json)
            # specimen_json["sample"].append(sample_json)
            # json_output["specimen"].append(specimen_json)

    with open(output, 'w', encoding='utf-8') as file:
        json.dump(json_output, file, ensure_ascii=False, indent=4)

    return readset_dict, sample_dict


def dna_raw_mean_coverage_check(sample, value, tumour):
    """ Mean Coverage DNA metric check """
    if not value:
        ret = "MISSING"
        logger.warning(f"Missing 'Mean Coverage' value for {sample}")
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
        logger.warning(f"Missing 'RNA Cluster' value for {sample}")
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
        logger.warning(f"Missing 'Dup. Rate (%)' value for {sample}")
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
        logger.warning(f"Missing 'Median Insert Size' value for {sample}")
    if float(value)<300:
        ret = "WARNING"
    elif float(value)<150:
        ret = "FAILED"
    else:
        ret = "PASS"
    return ret

def compute_md5(file_path, chunk_size=8 * 1024 * 1024):  # 8MB chunks
    """Compute or retrieve MD5 checksum of a file using EAFP style."""
    md5_file_path = f"{file_path}.md5"

    try:
        with open(md5_file_path, 'r') as f:
            line = f.readline()
            return line.split()[0]
    except (FileNotFoundError, IOError):
        pass  # Proceed to compute MD5 if .md5 file doesn't exist or can't be read

    # Compute MD5
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        while chunk := f.read(chunk_size):
            md5.update(chunk)
    return md5.hexdigest()

def main():
    """ Main """
    args = parse_arguments()

    if not args.output:
        output = f"{os.path.basename(args.input).split('.')[0]}.json"
    else:
        output = args.output
    if args.lane:
        lanes = list(args.lane)
    else:
        lanes = ["1", "2", "3", "4", "5", "6", "7", "8"]

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

if __name__ == '__main__':
    main()
