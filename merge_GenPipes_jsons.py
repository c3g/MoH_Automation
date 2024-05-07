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
    parser = argparse.ArgumentParser(prog=os.path.basename(__file__), description="Merges json files to be used as input for transferring GenPipes analysis.")
    parser.add_argument('-i', '--input', required=True, help="Json files from GenPipes analysis to be merged", nargs='+', type=open)
    parser.add_argument('-o', '--output', required=True, help="Output json filename.")
    args = parser.parse_args()

    merge_jsons(args.input, args.output)


def merge_jsons(json_files, output_file):
    """ Merges json files to be used as input for transferring GenPipes analysis. """
    merged_json = {
        "sample": []
    }
    merged_json_jobs = {}
    for json_file in json_files:
        data = json.load(json_file)
        samples = []
        for sample in data["sample"]:
            readsets = []
            for readset in sample["readset"]:
                if readset["readset_name"] not in merged_json_jobs:
                    merged_json_jobs[readset["readset_name"]] = []
                jobs = []
                for job in readset["job"]:
                    if job["job_status"] == "COMPLETED" and job["job_name"] not in merged_json_jobs[readset["readset_name"]]:
                        jobs.append(
                            {
                                "job_name": job["job_name"],
                                "file": job["file"]
                            }
                        )
                        merged_json_jobs[readset["readset_name"]].append(job["job_name"])
                if jobs:
                    readsets.append(
                        {
                            "readset_name": readset["readset_name"],
                            "job": jobs
                        }
                    )
            if readsets:
                samples.append(
                    {
                        "sample_name": sample["sample_name"],
                        "readset": readsets
                    }
                )
        for sample in samples:
            merged_json["sample"].append(sample)

    with open(output_file, 'w') as file:
        json.dump(merged_json, file, indent=4)


if __name__ == '__main__':
    main()
