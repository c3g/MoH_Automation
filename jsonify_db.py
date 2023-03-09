#!/usr/bin/env python3

import json
import glob
import sys
import re
import csv
import os
from datetime import datetime

def main():
    patient_dict = {}

    for sample in os.listdir("/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/raw_reads"):
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

    for run_file in glob.glob("/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/*.csv"):
        run_list = []
        with open(run_file, 'rt') as run_file_in:
            reader = csv.DictReader(run_file_in)
            for row in reader:
                run_list.append(row)
        json_output = {
                "bundle_config_uri": f"abacus:///lb/robot/research/processing/novaseq/20{run_list[0]['Processing Folder Name'][0:2]}/{run_list[0]['Processing Folder Name']}",
                "file_config_content": "",
                "file_config_type": "event",
                "project_fms_id": "",
                "project_name": "MOH-Q",
                "run_fms_id": "",
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
                    for _, element in enumerate(json_output["patient"]):
                        if f"{patient_dict[patient][sample][0]}" == element["patient_name"]:
                            existing_patient = element
                    # existing_patient = [element for _, element in enumerate(json_output["patient"]) if f"{patient_dict[patient][sample][0]}" == element["patient_name"]]
                    if existing_patient:
                        patient_json = existing_patient
                        # print(existing_patient)
                    else:
                        patient_json = {
                            "patient_fms_id": "",
                            "patient_name": f"{patient_dict[patient][sample][0]}",
                            "patient_cohort": f"{patient_dict[patient][sample][3]}",
                            "patient_institution": f"{patient_dict[patient][sample][4]}",
                            "sample": []
                            }

                    # print(patient_json)
                    existing_sample = None
                    for _, element in enumerate(patient_json["sample"]):
                        if f"{patient_dict[patient][sample][1]}" == element["sample_name"]:
                            existing_sample = element
                    # existing_sample = [element for _, element in enumerate(patient_json["sample"]) if f"{patient_dict[patient][sample][1]}" == element["sample_name"]]
                    if existing_sample:
                        sample_json = existing_sample
                    else:
                        if patient_dict[patient][sample][2].endswith("T"):
                            sample_tumour = True
                        else:
                            sample_tumour = False
                        sample_json = {
                            "sample_fms_id": "",
                            "sample_name": f"{patient_dict[patient][sample][1]}",
                            "sample_tumour": sample_tumour,
                            "readset": []
                            }

                    # transfer_folder = '/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*'
                    transfer_folder = 'transfer/*'
                    fastq1 = fastq2 = bundle_uri = ""
                    for filename in glob.glob(transfer_folder):
                        with open(filename, 'r') as file:
                            for line in file:
                                if patient_dict[patient][sample][1] in line:
                                    fields = line.split(",")
                                    bundle_uri = os.path.dirname(fields[0])
                                    if ".bam" in line:
                                        file_json = [
                                            {
                                                "file_content": f"{os.path.basename(fields[0])}"
                                                }
                                            ]
                                        break
                                    elif ".fastq" in line:
                                        if "R1.fastq" in line:
                                            fastq1 = os.path.basename(fields[0])
                                        elif "R2.fastq" in line:
                                            fastq2 = os.path.basename(fields[0])
                                        if fastq1 and fastq2:
                                            file_json = [
                                            {
                                                "file_content": f"{fastq1}",
                                                "file_extra_metadata": "R1"
                                                },
                                            {
                                                "file_content": f"{fastq2}",
                                                "file_extra_metadata": "R2"
                                                }
                                            ]
                                            break
                    metric_json = [
                        {
                            "metric_name": "raw_reads_count",
                            "metric_value": f"{run_row['Clusters']}"
                            },
                        {
                            "metric_name": "raw_duplication_rate",
                            "metric_value": f"{run_row['Dup. Rate (%)']}"
                            },
                        {
                            "metric_name": "raw_median_insert_size",
                            "metric_value": f"{run_row['Mapped Insert Size (median)']}"
                            },
                        {
                            "metric_name": "raw_mean_insert_size",
                            "metric_value": f"{run_row['Mapped Insert Size (mean)']}"
                            },
                        {
                            "metric_name": "raw_mean_coverage",
                            "metric_value": f"{run_row['Mean Coverage']}"
                            }
                        ]
                    # patient_json = {
                    #     "patient_fms_id": "",
                    #     "patient_name": f"{patient_dict[patient][sample][0]}",
                    #     "patient_cohort": f"{patient_dict[patient][sample][3]}",
                    #     "patient_institution": f"{patient_dict[patient][sample][4]}",
                    #     "sample": []
                    #     }
                    # if patient_dict[patient][sample][2].endswith("T"):
                    #     sample_tumour = True
                    # else:
                    #     sample_tumour = False
                    # sample_json = {
                    #     "sample_fms_id": "",
                    #     "sample_name": f"{patient_dict[patient][sample][1]}",
                    #     "sample_tumour": sample_tumour,
                    #     "readset": []
                    #     }
                    readset_json = {
                        "experiment_sequencing_technology": "",
                        "experiment_type": f"{run_row['Library Type']}",
                        "experiment_library_kit": "",
                        "experiment_kit_expiration_date": "",
                        "readset_name": f"{run_row['Sample ID']}.{run_row['Run ID']}_{run_row['Lane']}",
                        "readset_lane": f"{run_row['Lane']}",
                        "readset_adapter1": f"{run_row['i7 Adapter Sequence']}",
                        "readset_adapter2": f"{run_row['i5 Adapter Sequence']}",
                        "readset_sequencing_type": f"{run_row['Run Type']}",
                        "readset_quality_offset": "33",
                        "bundle_uri": bundle_uri,
                        "file": file_json,
                        "metric": metric_json
                        }
                    sample_json["readset"].append(readset_json)
                    patient_json["sample"].append(sample_json)
                    json_output["patient"].append(patient_json)
                except KeyError:
                    pass
                # json_output["patient"].append(
                #     {
                #         "patient_fms_id": "",
                #         "patient_name": f"{patient_dict[patient][sample][0]}",
                #         "patient_cohort": f"{patient_dict[patient][sample][3]}",
                #         "patient_institution": f"{patient_dict[patient][sample][4]}",
                #         "sample": [
                #             {
                #                 "sample_fms_id": "",
                #                 "sample_name": f"{patient_dict[patient][sample][1]}",
                #                 "sample_tumour": True,
                #                 "readset": [
                #                     {
                #                         "experiment_sequencing_technology": "",
                #                         "experiment_type": f"{run_row['Library Type']}",
                #                         "experiment_library_kit": "",
                #                         "experiment_kit_expiration_date": "",
                #                         "readset_name": f"{run_row['Sample ID']}.{run_row['Run ID']}_{run_row['Lane']}",
                #                         "readset_lane": f"{run_row['Lane']}",
                #                         "readset_adapter1": f"{run_row['i7 Adapter Sequence']}",
                #                         "readset_adapter2": f"{run_row['i5 Adapter Sequence']}",
                #                         "readset_sequencing_type": f"{run_row['Run Type']}",
                #                         "readset_quality_offset": "33",
                #                         "bundle_uri": bundle_uri,
                #                         "file": file_json,
                #                         "metric": [
                #                             {
                #                                 "metric_name": "raw_reads_count",
                #                                 "metric_value": f"{run_row['Clusters']}"
                #                             },
                #                             {
                #                                 "metric_name": "raw_duplication_rate",
                #                                 "metric_value": f"{run_row['Dup. Rate (%)']}"
                #                             },
                #                             {
                #                                 "metric_name": "raw_median_insert_size",
                #                                 "metric_value": f"{run_row['Mapped Insert Size (median)']}"
                #                             },
                #                             {
                #                                 "metric_name": "raw_mean_insert_size",
                #                                 "metric_value": f"{run_row['Mapped Insert Size (mean)']}"
                #                             },
                #                             {
                #                                 "metric_name": "raw_mean_coverage",
                #                                 "metric_value": f"{run_row['Mean Coverage']}"
                #                             }
                #                         ]
                #                     }
                #                 ]
                #             }
                #             ]
                #         }
                #     )
        # print(json.dumps(json_output, indent=4))
        with open(f"jsons/{run_list[0]['Processing Folder Name']}.json", 'w', encoding='utf-8') as f:
            json.dump(json_output, f, ensure_ascii=False, indent=4)

if __name__ == '__main__':
    main()
