#!/usr/bin/env python3

import argparse
import glob
import json
import os
import logging

logging.basicConfig(format='%(levelname)s: %(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


# Map the known prefixes to their path index offsets
PREFIX_OFFSETS = {
    "/lb/robot/research/processing/novaseq/20": 0,
    "/lb/robot/research/freezeman-processing/novaseqx/20": 0,
    "/lb/project/mugqic/projects/MOH/GQ_STAGING/mgc_mock_runs": 1,  # +1 indices for this prefix
}


def main():
    """ Main """
    parser = argparse.ArgumentParser(prog='transfer2json.py', description="Creates json file for project tracking database for a given transfer of data.")
    parser.add_argument('-i', '--input', required=True, help="Batch file from Globus.")
    parser.add_argument('-s', '--source', required=True, help="Source cluster of the transfer.")
    parser.add_argument('-d', '--destination', required=True, help="Cluster of destination for the transfer.")
    parser.add_argument('-o', '--output', required=False, help="Output json filename (Default: <input_filename>.json).")
    parser.add_argument('-j', '--genpipes', required=False, help="GenPipes json file when creating json for a GenPipes transfer.")
    parser.add_argument('--start', required=False, help="Start time of operation (format: YYYY-MM-DDTHH.MM.SS).")
    parser.add_argument('--stop', required=False, help="End time of operation (format: YYYY-MM-DDTHH.MM.SS).")
    parser.add_argument('--operation_cmd_line', required=True, help="Command used for transfer.")
    parser.add_argument('--delivery', required=False, help="Delivery json file.")
    args = parser.parse_args()

    if not args.output:
        base_name, _ = os.path.splitext(args.input)
        output = os.path.basename(base_name) + ".json"
    else:
        output = args.output

    if args.genpipes:
        jsonify_genpipes_transfer(
            batch_file = args.input,
            source = args.source,
            destination = args.destination.lower(),
            genpipes_json = args.genpipes,
            output = output,
            operation_cmd_line = args.operation_cmd_line,
            start = args.start,
            stop = args.stop
            )
    elif args.delivery:
        jsonify_delivery_transfer(
            batch_file = args.input,
            source = args.source,
            destination = args.destination.lower(),
            delivery_json = args.delivery,
            output = output,
            operation_cmd_line = args.operation_cmd_line,
            start = args.start,
            stop = args.stop
            )
    else:
        jsonify_run_processing_transfer(
            batch_file = args.input,
            source = args.source,
            destination = args.destination.lower(),
            output = output,
            operation_cmd_line = args.operation_cmd_line,
            start = args.start,
            stop = args.stop
            )

def jsonify_delivery_transfer(batch_file, source, destination, delivery_json, output, operation_cmd_line, start=None, stop=None):
    """Writing transfer json based on json delivery file"""
    with open(delivery_json, 'r') as json_file:
        delivery_json = json.load(json_file)
    delivery_file = {}
    for specimen in delivery_json['specimen']:
        for sample in specimen['sample']:
            for readset in sample['readset']:
                for file in readset['file']:
                    if file['name'] in delivery_file:
                        delivery_file[file['name']].append(readset['name'])
                    else:
                        delivery_file[file['name']] = [readset['name']]

    start = start.replace('.', ':').replace('T', ' ') if start else None
    stop = stop.replace('.', ':').replace('T', ' ') if stop else None
    json_output = {
        "operation_platform": source,
        "operation_cmd_line": operation_cmd_line,
        "job_start": start,
        "job_stop": stop,
        "readset": []
        }
    with open(batch_file, 'r') as file:
        for line in file:
            fields = line.split(" ")
            src_location_uri = f"{source}://{fields[0]}"
            dest_location_uri = f"{destination}://{fields[1].strip()}"
            current_file = os.path.basename(fields[0])
            if current_file in delivery_file:
                for readset_name in delivery_file[os.path.basename(current_file)]:
                    # relative_file_path = current_file.replace(fields[1], '')
                    # src_location_uri_file = f"{src_location_uri}{relative_file_path}"
                    # dest_location_uri_file = f"{dest_location_uri}{relative_file_path}"
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

    with open(output, 'w', encoding='utf-8') as file:
        json.dump(json_output, file, ensure_ascii=False, indent=4)


def _matched_prefix(line: str):
    """Return the matched prefix or None if not matched."""
    for pfx in PREFIX_OFFSETS.keys():
        if line.startswith(pfx):
            return pfx
    return None

def _append_file(json_output, readset_name, src_location_uri, dest_location_uri):
    """Append a file entry to an existing readset or create a new one."""
    for readset in json_output["readset"]:
        if readset["readset_name"] == readset_name:
            readset["file"].append(
                {
                    "src_location_uri": src_location_uri,
                    "dest_location_uri": dest_location_uri
                }
            )
            return
    # Not found: create new
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

def jsonify_run_processing_transfer(batch_file, source, destination, output, operation_cmd_line, start=None, stop=None):
    """Writing transfer json based on batch file with path-index offsets per prefix."""
    start = start.replace('.', ':').replace('T', ' ') if start else None
    stop = stop.replace('.', ':').replace('T', ' ') if stop else None

    json_output = {
        "operation_platform": source,
        "operation_cmd_line": operation_cmd_line,
        "job_start": start,
        "job_stop": stop,
        "readset": []
    }

    with open(batch_file, 'r') as file:
        for raw_line in file:
            line = raw_line.strip()
            pfx = _matched_prefix(line)
            if not pfx:
                continue # skip lines that are not one of the known prefixes

            fields = line.split()
            if len(fields) < 2:
                continue # malformed line

            src_path = fields[0]
            dest_path = fields[1]

            path_l = src_path.split('/')
            offset = PREFIX_OFFSETS[pfx]

            # Build URIs once
            src_location_uri = f"abacus://{src_path}"
            dest_location_uri = f"{destination}://{dest_path.strip()}"

            # BAM CASE
            if ".bam" in src_path:
                try:
                    readset_name = f"{path_l[10 + offset]}.{path_l[11 + offset].replace('run', '')}"
                except IndexError:
                    continue

                _append_file(json_output, readset_name, src_location_uri, dest_location_uri)

            # FASTQ CASE
            elif ".fastq" in src_path:
                try:
                    run_part = path_l[7 + offset]
                    lane = path_l[8 + offset].split('.')[1]
                    sample_name = path_l[10 + offset].split('_')[1]
                except (IndexError, ValueError):
                    continue

                is_mgc = pfx.endswith("mgc_mock_runs")

                if is_mgc:
                    try:
                        _, instrument_full = run_part.split('_')
                        instrument = instrument_full.split('-')[0]
                    except ValueError:
                        continue

                    readset_name = f"{sample_name}.{instrument}_{lane}"

                else:
                    # regular case
                    parts = run_part.split('_')
                    if len(parts) < 3:
                        continue
                    instrument = parts[1]
                    runnum = parts[2]
                    readset_name = f"{sample_name}.{instrument}_{runnum}_{lane}"

                _append_file(json_output, readset_name, src_location_uri, dest_location_uri)

    with open(output, 'w', encoding='utf-8') as file:
        json.dump(json_output, file, ensure_ascii=False, indent=4)


def jsonify_genpipes_transfer(batch_file, source, destination, genpipes_json, output, operation_cmd_line, start=None, stop=None):
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
    start = start.replace('.', ':').replace('T', ' ') if start else None
    stop = stop.replace('.', ':').replace('T', ' ') if stop else None
    json_output = {
        "operation_platform": source,
        "operation_cmd_line": operation_cmd_line,
        "job_start": start,
        "job_stop": stop,
        "readset": []
        }
    with open(batch_file, 'r') as file:
        for line in file:
            fields = line.split(" ")
            if line.startswith("--recursive"):
                src_location_uri = f"{source}://{fields[1]}"
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
                src_location_uri = f"{source}://{fields[0]}"
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

    with open(output, 'w', encoding='utf-8') as file:
        json.dump(json_output, file, ensure_ascii=False, indent=4)

if __name__ == '__main__':
    main()
