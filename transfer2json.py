#!/usr/bin/env python3

import argparse
import json
import os
import logging

logging.basicConfig(format='%(levelname)s: %(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    """ Main """
    parser = argparse.ArgumentParser(prog='transfer2json.py', description="Creates json file for project tracking database for a given transfer of data.")
    parser.add_argument('-i', '--input', required=True, help="Input align_bwa_mem.csv file from Run Processing.")
    parser.add_argument('-o', '--output', required=False, help="Output json filename (Default: <input_filename>.json).")
    parser.add_argument('--operation_cmd_line', required=True, help="Command used for transfer.")
    args = parser.parse_args()

    if not args.output:
        output = os.path.basename(args.input).replace('.txt', '.json')
    else:
        output = args.output

    # jsonify_run_processing_transfer(args.input, output)
    jsonify_run_processing_transfer(args.input, output, args.operation_cmd_line)

def jsonify_run_processing_transfer(batch_file, output, operation_cmd_line):
    """Writing transfer json based on batch file"""
    transfer_dict = {}
    json_output = {
        "operation_platform": "abacus",
        "operation_cmd_line": operation_cmd_line,
        "readset": []
        }
    with open(batch_file, 'r') as file:
        for line in file:
            if line.startswith("/lb/robot/research/processing/novaseq/20"):
                fields = line.split(" ")
                if ".bam" in line:
                    readset_name = f"{fields[0].split('/')[10]}_{fields[0].split('/')[11]}"
                    src_location_uri = f"abacus://{fields[0]}"
                    dest_location_uri = f"beluga://{fields[1].strip()}"
                    if "bam.bai" in line:
                        if readset_name in transfer_dict:
                            transfer_dict[readset_name]["bai_dest_location_uri"] = dest_location_uri
                        else:
                            transfer_dict[readset_name] = {
                                "bai_dest_location_uri": dest_location_uri
                            }
                        transfer_dict[readset_name]["bai_src_location_uri"] = src_location_uri
                    else:
                        if readset_name in transfer_dict:
                            transfer_dict[readset_name]["bam_dest_location_uri"] = dest_location_uri
                        else:
                            transfer_dict[readset_name] = {
                                "bam_dest_location_uri": dest_location_uri
                            }
                        transfer_dict[readset_name]["bam_src_location_uri"] = src_location_uri

                elif ".fastq" in line:
                    readset_name = fields[0].split('/')[10].replace("Sample_", "new")
                    src_location_uri = f"abacus://{fields[0]}"
                    dest_location_uri = f"beluga://{fields[1].strip()}"
                    if "_R1_" in line:
                        if readset_name in transfer_dict:
                            transfer_dict[readset_name]["fastq1_dest_location_uri"] = dest_location_uri
                        else:
                            transfer_dict[readset_name] = {
                                "fastq1_dest_location_uri": dest_location_uri
                            }
                        transfer_dict[readset_name]["fastq1_src_location_uri"] = src_location_uri
                    elif "_R2_" in line:
                        if readset_name in transfer_dict:
                            transfer_dict[readset_name]["fastq2_dest_location_uri"] = dest_location_uri
                        else:
                            transfer_dict[readset_name] = {
                                "fastq2_dest_location_uri": dest_location_uri
                            }
                        transfer_dict[readset_name]["fastq2_src_location_uri"] = src_location_uri

    for readset in transfer_dict:
        if "bam_src_location_uri" in transfer_dict[readset]:
            file_json = [
                {
                    "src_location_uri": transfer_dict[readset]["bam_src_location_uri"],
                    "dest_location_uri": transfer_dict[readset]["bam_dest_location_uri"]
                },
                {
                    "src_location_uri": transfer_dict[readset]["bai_src_location_uri"],
                    "dest_location_uri": transfer_dict[readset]["bai_dest_location_uri"],

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
    with open(output, 'w', encoding='utf-8') as file:
        json.dump(json_output, file, ensure_ascii=False, indent=4)


if __name__ == '__main__':
    main()
