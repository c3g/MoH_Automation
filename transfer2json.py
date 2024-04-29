#!/usr/bin/env python3

import argparse
import glob
import json
import os
import logging

logging.basicConfig(format='%(levelname)s: %(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    """ Main """
    parser = argparse.ArgumentParser(prog='transfer2json.py', description="Creates json file for project tracking database for a given transfer of data.")
    parser.add_argument('-i', '--input', required=True, help="Input align_bwa_mem.csv file from Run Processing.")
    parser.add_argument('-d', '--destination', required=True, help="Cluster of destination for the transfer.")
    parser.add_argument('-o', '--output', required=False, help="Output json filename (Default: <input_filename>.json).")
    parser.add_argument('-j', '--genpipes', required=False, help="GenPipes json file when creating json for a GenPipes transfer.")
    parser.add_argument('--operation_cmd_line', required=True, help="Command used for transfer.")
    args = parser.parse_args()

    if not args.output:
        base_name, _ = os.path.splitext(args.input)
        output = os.path.basename(base_name) + ".json"
    else:
        output = args.output

    if args.genpipes:
        jsonify_genpipes_transfer(args.input, args.destination.lower(), args.genpipes, output, args.operation_cmd_line)
    else:
        jsonify_run_processing_transfer(args.input, args.destination.lower(), output, args.operation_cmd_line)

def jsonify_run_processing_transfer(batch_file, destination, output, operation_cmd_line):
    """Writing transfer json based on batch file"""
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
                    path_l = fields[0].split('/')
                    readset_name = f"{path_l[10]}.{path_l[11].replace('run', '')}"
                    src_location_uri = f"abacus://{fields[0]}"
                    dest_location_uri = f"{destination}://{fields[1].strip()}"
                    if readset_name in [readset["readset_name"] for readset in json_output["readset"]]:
                        for readset in json_output["readset"]:
                            if readset_name == readset["readset_name"]:
                                readset["file"].append(
                                    {
                                        "src_location_uri": src_location_uri,
                                        "dest_location_uri": dest_location_uri
                                    }
                                )
                    else:
                        json_output["readset"].append(
                            {
                                "readset_name": readset_name,
                                "file": [
                                    {
                                        "src_location_uri": src_location_uri,
                                        "dest_location_uri": dest_location_uri
                                    }
                                ]
                            }
                        )

                elif ".fastq" in line:
                    path_l = fields[0].split('/')
                    run_name_l = path_l[7].split('_')
                    sample_name = path_l[10].split('_')[1]
                    lane = path_l[8].split('.')[1]
                    readset_name = f"{sample_name}.{run_name_l[1]}_{run_name_l[2]}_{lane}"
                    src_location_uri = f"abacus://{fields[0]}"
                    dest_location_uri = f"{destination}://{fields[1].strip()}"
                    if readset_name in [readset["readset_name"] for readset in json_output["readset"]]:
                        for readset in json_output["readset"]:
                            if readset_name == readset["readset_name"]:
                                readset["file"].append(
                                    {
                                        "src_location_uri": src_location_uri,
                                        "dest_location_uri": dest_location_uri
                                    }
                                )
                    else:
                        json_output["readset"].append(
                            {
                                "readset_name": readset_name,
                                "file": [
                                    {
                                        "src_location_uri": src_location_uri,
                                        "dest_location_uri": dest_location_uri
                                    }
                                ]
                            }
                        )

    # print(json.dumps(json_output, ensure_ascii=False, indent=4))

    with open(output, 'w', encoding='utf-8') as file:
        json.dump(json_output, file, ensure_ascii=False, indent=4)


def jsonify_genpipes_transfer(batch_file, destination, genpipes_json, output, operation_cmd_line):
    """Writing transfer json based on batch file"""
    with open(genpipes_json, 'r') as json_file:
        genpipes_json = json.load(json_file)
    genpipes_file = {}
    for sample in genpipes_json['sample']:
        for readset in sample['readset']:
            for job in readset['job']:
                for file in job['file']:
                    if file['file_name'] in genpipes_file:
                        genpipes_file[file['file_name']].append(readset['readset_name'])
                    else:
                        genpipes_file[file['file_name']] = [readset['readset_name']]
    json_output = {
        "operation_platform": "abacus",
        "operation_cmd_line": operation_cmd_line,
        "readset": []
        }
    with open(batch_file, 'r') as file:
        for line in file:
            fields = line.split(" ")
            if line.startswith("--recursive"):
                src_location_uri = f"abacus://{fields[1]}"
                dest_location_uri = f"{destination}://{fields[2].strip()}"
                filename = glob.glob(os.path.join(fields[1], '**'), recursive=True)
                for current_file in filename:
                    if os.path.basename(current_file) in genpipes_file:
                        for readset_name in genpipes_file[os.path.basename(current_file)]:
                            relative_file_path = current_file.replace(fields[1], '')
                            src_location_uri_file = f"{src_location_uri}{relative_file_path}"
                            dest_location_uri_file = f"{dest_location_uri}{relative_file_path}"
                            if readset_name in [readset["readset_name"] for readset in json_output["readset"]]:
                                for readset in json_output["readset"]:
                                    if readset_name == readset["readset_name"]:
                                        readset["file"].append(
                                            {
                                                "src_location_uri": src_location_uri_file,
                                                "dest_location_uri": dest_location_uri_file
                                            }
                                        )
                            else:
                                json_output["readset"].append(
                                    {
                                        "readset_name": readset_name,
                                        "file": [
                                            {
                                                "src_location_uri": src_location_uri_file,
                                                "dest_location_uri": dest_location_uri_file
                                            }
                                        ]
                                    }
                                )
            else:
                src_location_uri = f"abacus://{fields[0]}"
                dest_location_uri = f"{destination}://{fields[1].strip()}"
                current_file = os.path.basename(fields[0])
                if current_file in genpipes_file:
                    for readset_name in genpipes_file[os.path.basename(current_file)]:
                        relative_file_path = current_file.replace(fields[1], '')
                        src_location_uri_file = f"{src_location_uri}{relative_file_path}"
                        dest_location_uri_file = f"{dest_location_uri}{relative_file_path}"
                        if readset_name in [readset["readset_name"] for readset in json_output["readset"]]:
                            for readset in json_output["readset"]:
                                if readset_name == readset["readset_name"]:
                                    readset["file"].append(
                                        {
                                            "src_location_uri": src_location_uri_file,
                                            "dest_location_uri": dest_location_uri_file
                                        }
                                    )
                        else:
                            json_output["readset"].append(
                                {
                                    "readset_name": readset_name,
                                    "file": [
                                        {
                                            "src_location_uri": src_location_uri_file,
                                            "dest_location_uri": dest_location_uri_file
                                        }
                                    ]
                                }
                            )

    # print(json.dumps(json_output, ensure_ascii=False, indent=4))

    with open(output, 'w', encoding='utf-8') as file:
        json.dump(json_output, file, ensure_ascii=False, indent=4)

if __name__ == '__main__':
    main()
