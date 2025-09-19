#!/usr/bin/env python3

# standard import
import argparse
import json
import datetime
import logging
import os
import re
import csv
from collections import defaultdict
from io import StringIO
from pathlib import Path
import signal
import sys
import time

import botocore
import boto3
import globus_sdk
import markdown
from pymdownx import emoji
from bs4 import BeautifulSoup

from dateutil.parser import parse
from retry import retry

extensions = [
    'markdown.extensions.tables',
    'pymdownx.magiclink',
    'pymdownx.betterem',
    'pymdownx.tilde',
    'pymdownx.emoji',
    'pymdownx.tasklist',
    'pymdownx.superfences',
    'pymdownx.saneheaders',
    'footnotes'
]

extension_configs = {
    "pymdownx.magiclink": {
        "repo_url_shortener": True,
        "repo_url_shorthand": True,
        "provider": "github",
        "user": "facelessuser",
        "repo": "pymdown-extensions"
    },
    "pymdownx.tilde": {
        "subscript": False
    },
    "pymdownx.emoji": {
        "emoji_index": emoji.emojione,
        "emoji_generator": emoji.to_png_sprite,
    }
}

logger = logging.getLogger(__name__)

def main():
    """
    Main function to prepare delivery to the bucket.
    """
    parser = argparse.ArgumentParser(prog=__name__, description="Prepare delivery to the bucket.")
    parser.add_argument('-i', '--input', required=True, help="The json file returned by pt_cli digest delivery.")
    parser.add_argument('-l', '--list_file', required=True, help="The name of the list file for transferred files.")
    # parser.add_argument('-f', '--force', required=False, help="Ignores thresholds over metrics to deliver.", action='store_true')
    parser.add_argument('--ignore_alignment', required=False, help="Ignores alignments for delivery; aka doesn't deliver it.", action='store_true')
    parser.add_argument('--update_metrics', required=False, help="Forces Key_metrics.csv and Warnings.html files generation.", action='store_true')
    parser.add_argument('--update_methods', required=False, help="Forces Methods.html file generation.", action='store_true')
    parser.add_argument('--update_readme', required=False, help="Forces Readme.html file generation.", action='store_true')
    parser.add_argument('--loglevel', help='Sets logging level', choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'], default='INFO')
    args = parser.parse_args()

    log_level = getattr(logging, args.loglevel.upper(), None)

    if log_level is not None:
        logger.setLevel(log_level)
        logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    else:
        logger.setLevel(logging.INFO)
        logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')

    # Path to your environment variables file
    globus_env_file_path = os.path.expanduser('~/.config/globus_cli/globus_cli_env.sh')
    boto3_env_file_path = os.path.expanduser('~/.config/openstack/boto3_env.sh')

    # Read the environment variables from the file
    read_env_file(globus_env_file_path)
    read_env_file(boto3_env_file_path)

    # Now you can access the environment variables
    client_id = os.getenv('GLOBUS_CLI_CLIENT_ID')
    client_secret = os.getenv('GLOBUS_CLI_CLIENT_SECRET')
    secret_id = os.getenv('OPENSTACK_SECRET_ID')
    secret_key = os.getenv('OPENSTACK_SECRET_KEY')

    # Initialize the Globus SDK client
    client = globus_sdk.ConfidentialAppAuthClient(client_id, client_secret)
    token_response = client.oauth2_client_credentials_tokens()

    # Initialize the S3 client
    s3_client = boto3.client('s3',
        endpoint_url='https://objets.juno.calculquebec.ca',
        aws_access_key_id=secret_id,
        aws_secret_access_key=secret_key,
        verify=True
        )

    # S3 bucket details
    bucket_name = 'MOH-Q'

    transfer_tokens = token_response.by_resource_server['transfer.api.globus.org']
    transfer_authorizer = globus_sdk.AccessTokenAuthorizer(transfer_tokens['access_token'])
    transfer_client = globus_sdk.TransferClient(authorizer=transfer_authorizer)


    with open(args.input, 'r') as json_file:
        json_content = json.load(json_file)

    experiment_nucleic_acid_type = json_content["experiment_nucleic_acid_type"]
    location_endpoint = json_content["location_endpoint"]
    operation_list = json_content["operation"]

    now = datetime.datetime.today().strftime("%Y-%m-%dT%H.%M.%S")

    for patient in json_content["specimen"]:
        patient_name = patient["name"]
        cohort = patient["cohort"]
        institution = patient["institution"]
        sample_name_dna_n = sample_name_dna_t = sample_name_rna = patient_name
        dna = rna = False
        for sample in patient["sample"]:
            sample_name = sample["name"]
            tumour = sample["tumour"]
            if experiment_nucleic_acid_type == "DNA":
                dna = True
                if tumour:
                    sample_name_dna_t = sample_name
                else:
                    sample_name_dna_n = sample_name
                dedup_coverage = extract_metrics(sample, "dedup_coverage")
                if not isinstance(dedup_coverage, list):
                    if float(dedup_coverage) <= 30 and not tumour:
                        logger.warning(f"Metric dedup_coverage for Sample {sample_name} is <= 30")
                        # if not args.force:
                            # dna = False
                    if float(dedup_coverage) <= 80 and tumour:
                        logger.warning(f"Metric dedup_coverage for Sample {sample_name} is <= 80")
                        # if not args.force:
                            # dna = False
                else:
                    logger.warning(f"Multiple dedup_coverage values for Sample {sample_name}: {dedup_coverage}.")
                    # if not args.force:
                    #     logger.warning(f"Skipping delivery for Sample {sample_name}. To force delivery, use '-f' option.")
                    #     dna = False

            elif experiment_nucleic_acid_type == "RNA":
                rna = True
                sample_name_rna = sample_name
                raw_reads_count = extract_metrics(sample, "raw_reads_count")
                if not isinstance(raw_reads_count, list):
                    if float(raw_reads_count) <= 80000000:
                        logger.warning(f"Metric raw_reads_count for Sample {sample_name} is <= 80.000.000")
                        # if not args.force:
                            # rna = False
                else:
                    logger.warning(f"Multiple raw_reads_count values for Sample {sample_name}: {raw_reads_count}.")
                    # if not args.force:
                    #     logger.warning(f"Skipping delivery for Sample {sample_name}. To force delivery, use '-f' option.")
                    #     rna = False

        # Folders used for Delivery
        with open(f"{os.path.join(os.path.dirname(os.path.abspath(__file__)), 'globus_collections.json')}", "r") as globus_collection_file:
            globus_collection = json.load(globus_collection_file)

        # Input UUID
        in_uuid = globus_collection["robot_endpoints"][location_endpoint]["uuid"]
        # Input Folder
        in_base_path = globus_collection["robot_endpoints"][location_endpoint]["base_path"]
        if location_endpoint == "abacus":
            # Abacus RawData is located somewhere else
            in_uuid_abacus_rawdata = globus_collection["robot_endpoints"]["abacus_rawdata"]["uuid"]
            in_base_path_abacus_rawdata = globus_collection["robot_endpoints"]["abacus_rawdata"]["base_path"]
        else:
            in_uuid_abacus_rawdata = None
            in_base_path_abacus_rawdata = None
        # Output UUID
        out_uuid = globus_collection["robot_endpoints"]["sd4h"]["uuid"]
        # Output Folder
        out_base_path = globus_collection["robot_endpoints"]["sd4h"]["base_path"]
        # Contains Warnings.txt Readme.txt Log.txt and all subfolders
        out_folder = os.path.join(out_base_path, institution, cohort, patient_name)
        # Contains raw bams and fastqs
        raw_folder = os.path.join(out_folder, "raw_data")
        # Contains all variants the subfolder
        var_folder = os.path.join(out_folder, "variants")
        # Contains all the vcfs from the callers
        cal_folder = os.path.join(var_folder, "caller_vcfs")
        # Contains raw cnv
        raw_cnv_folder = os.path.join(out_folder, "raw_cnv")
        # Contains the analysis bams
        alignment_folder = os.path.join(out_folder, "alignment")
        # Contains the ini files
        param_folder = os.path.join(out_folder, "parameters")
        # Contains the reports
        reports_folder = os.path.join(out_folder, "reports")
        pcgr_folder = os.path.join(reports_folder, "pcgr")
        # Contains expression from Kallisto
        expression_folder = os.path.join(out_folder, "expression")
        # Contains all svariants
        svar_folder = os.path.join(out_folder, "svariants")
        linx_folder = os.path.join(svar_folder, "linx")

        warning_file = os.path.join(out_folder, "Warnings.html")
        readme_file = os.path.join(out_folder, "Readme.html")
        methods_file = os.path.join(out_folder, "Methods.html")
        key_metrics_file = os.path.join(out_folder, "Key_metrics.csv")
        file_dict = {}
        # For Abacus rawdata special handling
        file_dict_rawdata = {}
        ini_dict = {}
        metrics_dict = defaultdict(lambda: None)
        # Check if the key metrics file exists
        key_metrics_file_outpath = remove_path_parts(key_metrics_file, out_base_path)
        if file_exists(s3_client, bucket_name, key_metrics_file_outpath):
            response = s3_client.get_object(Bucket=bucket_name, Key=key_metrics_file_outpath)
            file_content = response['Body'].read().decode('utf-8')
            existing_metrics = parse_csv(file_content)
            if existing_metrics and 'Readset' not in existing_metrics[0]:
                for row in existing_metrics:
                    row['Readset'] = 'NA'
        else:
            existing_metrics = []
        # Check if the warning file exists
        warning_file_outpath = remove_path_parts(warning_file, out_base_path)
        if file_exists(s3_client, bucket_name, warning_file_outpath):
            response = s3_client.get_object(Bucket=bucket_name, Key=warning_file_outpath)
            warnings_content = response['Body'].read().decode('utf-8')
        else:
            warnings_content = ""

        # Populates dna data
        if dna:
            deliver_dna(
                args.ignore_alignment,
                raw_folder,
                var_folder,
                svar_folder,
                linx_folder,
                cal_folder,
                raw_cnv_folder,
                alignment_folder,
                reports_folder,
                pcgr_folder,
                patient,
                out_base_path,
                in_base_path,
                in_base_path_abacus_rawdata,
                location_endpoint,
                file_dict,
                file_dict_rawdata,
                metrics_dict
                )
            for index, operation_object in enumerate(operation_list):
                ini_dict[os.path.join(param_folder, f"{patient_name}.TumourPair.{index+1}.ini")] = operation_object["config_data"]


        # Populates rna data
        elif rna:
            deliver_rna(
                args.ignore_alignment,
                raw_folder,
                expression_folder,
                var_folder,
                alignment_folder,
                reports_folder,
                pcgr_folder,
                patient,
                out_base_path,
                in_base_path,
                in_base_path_abacus_rawdata,
                location_endpoint,
                file_dict,
                file_dict_rawdata,
                metrics_dict
                )
            for index, operation_object in enumerate(operation_list):
                match = re.search(r'-t (\w+)', operation_object["cmd_line"])
                protocol = match.group(1) if match else None
                if protocol == "cancer":
                    ini_dict[os.path.join(param_folder, f"{patient_name}.RNA.Variants.{index+1}.ini")] = operation_object["config_data"]
                else:
                    ini_dict[os.path.join(param_folder, f"{patient_name}.RNA.Light.{index+1}.ini")] = operation_object["config_data"]

        # Not implemented yet
        # if rna and dna:
        #     os.makedirs(var_folder, exist_ok=True)
        #     final_vcf = extract_fileloc_field(connection, sample.sample, "Final_VCF")
        #     updated = get_link_log(final_vcf, var_folder, f"{sample.sample_true}.vcf.gz", log, updated, old_log)

        already_delivered_files = list_s3_files(s3_client, bucket_name, remove_path_parts(out_folder, out_base_path))
        # for src_file, dest_file in file_dict.items():
        #     print(f"Transferring {src_file} to {dest_file}")

        transfer_label = f"Delivery_{patient_name}_{experiment_nucleic_acid_type}_{now}"
        logger.debug(f"Transfer label: {transfer_label}")

        # list_file = os.path.join(args.list_file, f"Delivery_{patient_name}_{experiment_nucleic_acid_type}_{now}.list")

        transferred_files = []

        # Transfer for Abacus rawdata only
        if file_dict_rawdata:
            logger.debug(f"Transferring Abacus rawdata files: {file_dict_rawdata}")
            task_id = transfer_files_with_sync(in_uuid_abacus_rawdata, out_uuid, file_dict_rawdata, transfer_client, f"{transfer_label}_rawdata")
            display_transfer_status(transfer_client, task_id, s3_client, bucket_name)
            transferred_files, _ = get_transfer_event_log(transfer_client, task_id, s3_client, bucket_name)
            already_delivered_files.extend(format_timestamps(transferred_files))

        # Transfer for all endpoints (including abacus after rawdata split)
        if file_dict:
            logger.debug(f"Transferring files: {file_dict}")
            task_id = transfer_files_with_sync(in_uuid, out_uuid, file_dict, transfer_client, transfer_label)
            display_transfer_status(transfer_client, task_id, s3_client, bucket_name)
            transferred_files, _ = get_transfer_event_log(transfer_client, task_id, s3_client, bucket_name)
            already_delivered_files.extend(format_timestamps(transferred_files))

        # task_id = transfer_files_with_sync(in_uuid, out_uuid, file_dict, transfer_client, transfer_label)
        # display_transfer_status(transfer_client, task_id, s3_client, bucket_name)
        # transferred_files, _ = get_transfer_event_log(transfer_client, task_id, s3_client, bucket_name)
        # transferred_files = format_timestamps(transferred_files)
        # already_delivered_files.extend(transferred_files)
        all_delivered_files = list(set(already_delivered_files))

        if transferred_files:
            lines = []
            for file, _ in transferred_files:
                for src_file, dest_file in file_dict.items():
                    if dest_file == file:
                        lines.append(f"{os.path.join(in_base_path, src_file)} {os.path.join(out_base_path, dest_file)}")
                for src_file, dest_file in file_dict_rawdata.items():
                    if dest_file == file:
                        lines.append(f"{os.path.join(in_base_path_abacus_rawdata, src_file)} {os.path.join(out_base_path, dest_file)}")
            with open(args.list_file, 'w') as list_file:
                list_file.write('\n'.join(lines))
            for ini_file_name, ini_content in ini_dict.items():
                s3_client.put_object(Bucket=bucket_name, Key=remove_path_parts(ini_file_name, out_base_path), Body=ini_content)
            # Read already delivered warnings
            existing_warnings_dict, soup = parse_existing_warnings(warnings_content)
            updated_warnings_dict = update_warnings(existing_warnings_dict, metrics_dict, soup)
            updated_warnings_content = generate_updated_warnings_content(updated_warnings_dict, soup)
            s3_client.put_object(Bucket=bucket_name, Key=remove_path_parts(warning_file, out_base_path), Body=updated_warnings_content)
        else:
            logger.warning("No new file were transferred.")

        # Add methods file
        if args.update_methods or transferred_files:
            methods_content = generate_methods()
            s3_client.put_object(Bucket=bucket_name, Key=remove_path_parts(methods_file, out_base_path), Body=methods_content)

        if args.update_metrics or transferred_files:
            # Read already delivered key metrics
            updated_metrics_list = update_metrics(existing_metrics, metrics_dict)
            metrics_content = generate_csv_content(updated_metrics_list)
            s3_client.put_object(Bucket=bucket_name, Key=remove_path_parts(key_metrics_file, out_base_path), Body=metrics_content)

        if args.update_readme or transferred_files:
            readme_content = generate_readme(patient_name, sample_name_dna_n, sample_name_dna_t, sample_name_rna, all_delivered_files)
            s3_client.put_object(Bucket=bucket_name, Key=remove_path_parts(readme_file, out_base_path), Body=readme_content)


def extract_metrics(sample_content, metric_name):
    """Extracts a metric from a sample content dictionary.
    
    Returns a float if all values are the same, otherwise returns a list of unique float values.
    """
    metric_values = [
        metric["value"]
        for readset in sample_content["readset"]
        for metric in readset["metric"]
        if metric["name"] == metric_name
    ]

    if not metric_values:
        raise ValueError(f"No values found for metric '{metric_name}'.")

    unique_values = sorted(set(float(v) for v in metric_values))

    if len(unique_values) == 1:
        return unique_values[0]

    return unique_values


def deliver_dna(
    ignore_alignment,
    raw_folder,
    var_folder,
    svar_folder,
    linx_folder,
    cal_folder,
    raw_cnv_folder,
    alignment_folder,
    reports_folder,
    pcgr_folder,
    patient,
    out_base_path,
    in_base_path,
    in_base_path_abacus_rawdata,
    location_endpoint,
    file_dict,
    file_dict_rawdata,
    metrics_dict
    ):

    metric_mapping = {
        "raw_mean_coverage": "Raw_Mean_Coverage",
        "raw_duplication_rate": "Raw_Duplication_Rate",
        "raw_median_insert_size": "Raw_Median_Insert_Size",
        "raw_mean_insert_size": "Raw_Mean_Insert_Size",
        "raw_reads_count": "Raw_Reads_Count",
        "median_insert_size": "Median_Insert_Size",
        "bases_over_q30_percent": "WGS_Bases_Over_Q30",
        "aligned_reads_count": "WGS_Min_Aligned_Reads_Delivered",
        "dedup_coverage": "WGS_Dedup_Coverage",
        "contamination": "WGS_Contamination",
        "concordance": "Concordance",
        "purity": "Purity"
    }

    # Compile the regex patterns
    variant_pattern = re.compile(r"\.ensemble\.(germline|somatic)\.vt\.annot\.vcf\.gz")
    cal_pattern = re.compile(r"(mutect2|strelka2|vardict|varscan2)\.(somatic|germline)\.vt\.vcf\.gz")
    svar_pattern = re.compile(
        r"(gridss|gripss\.filtered\.(somatic|germline))\.vcf\.gz|driver\.catalog\.(somatic|germline)\.tsv|circos\.png$|\.purple_(ensemble|sv)\.zip$"
    )
    reports_pattern = re.compile(r"(\.multiqc|\.pcgr_acmg\.grch38\.flexdb)\.html$|\.(cpsr|pcgr)\.zip$")
    pcgr_pattern = re.compile(r"(\.pcgr_acmg\.grch38\.(maf|snvs_indels\.tiers\.tsv|cna_segments\.tsv\.gz))$")

    for sample in patient["sample"]:
        sample_name = sample["name"]
        for readset in sample["readset"]:
            readset_name = readset["name"]
            for file in readset["file"]:
                file_name = file["name"]
                # Handle Abacus rawdata path located in /lb/robot/research/freezeman-processing/novaseqx/
                if location_endpoint == "abacus" and file["location"].startswith(in_base_path_abacus_rawdata):
                    file_location = remove_path_parts(file["location"], in_base_path_abacus_rawdata)
                    # To workaround issue with RP naming regarding raw data being non unique we have to use file_location for file_name
                    file_dict_rawdata[file_location] = os.path.join(remove_path_parts(raw_folder, out_base_path), file_name)
                    continue
                # Usual case
                file_location = remove_path_parts(file["location"], in_base_path)
                # raw_data
                if "MAIN/raw_reads/" in file_location or "GQ_STAGING" in file_location:
                    # To workaround issue with RP naming regarding raw data being non unique we have to use file_location for file_name
                    file_dict[file_location] = os.path.join(remove_path_parts(raw_folder, out_base_path), os.path.basename(file_location))
                # variants
                elif variant_pattern.search(file_name):
                    file_dict[file_location] = os.path.join(remove_path_parts(var_folder, out_base_path), file_name)
                # variants/caller_vcfs
                elif cal_pattern.search(file_name):
                    file_dict[file_location] = os.path.join(remove_path_parts(cal_folder, out_base_path), file_name)
                # raw_cnv
                elif "cnvkit.vcf.gz" in file_name:
                    file_dict[file_location] = os.path.join(remove_path_parts(raw_cnv_folder, out_base_path), file_name)
                # alignment
                elif "MAIN/alignment/" in file_location and not ignore_alignment:
                    # Rename to remove extra part from dna bams
                    file_name = file_name.replace(".sorted.dup.recal", "")
                    file_dict[file_location] = os.path.join(remove_path_parts(alignment_folder, out_base_path), file_name)
                # svariants
                elif svar_pattern.search(file_name):
                    file_dict[file_location] = os.path.join(remove_path_parts(svar_folder, out_base_path), file_name)
                # svariants/linx
                elif "linx" in file_location:
                    file_dict[file_location] = os.path.join(remove_path_parts(linx_folder, out_base_path), file_name)
                # reports
                elif reports_pattern.search(file_name):
                    match = reports_pattern.search(file_name)
                    if match:
                        # Insert '_D' before the matched part to differentiate between DNA and RNA reports
                        start = match.start()
                        file_name = file_name[:start] + "_D" + file_name[start:]
                    # Rename pcgr report
                    if "pcgr_acmg" in file_name:
                        file_name = file_name.replace("_acmg.grch38.flexdb", "")
                    file_dict[file_location] = os.path.join(remove_path_parts(reports_folder, out_base_path), file_name)
                # reports/pcgr
                elif pcgr_pattern.search(file_name):
                    match = pcgr_pattern.search(file_name)
                    if match:
                        # Insert '_D' before the matched part to differentiate between DNA and RNA reports
                        start = match.start()
                        file_name = file_name[:start] + "_D" + file_name[start:]
                    file_dict[file_location] = os.path.join(remove_path_parts(pcgr_folder, out_base_path), file_name)
            for metric in readset["metric"]:
                metric_name = metric["name"]
                metric_value = metric["value"]
                metric_flag = metric["flag"]
                # metric_aggregate = metric["aggregate"]
                # update_metrics(metrics_dict, metric_mapping, sample_name, readset_name, metric_name, metric_value, metric_flag)
                update_metrics_dict(metrics_dict, metric_mapping, sample_name, readset_name, metric_name, metric_value, metric_flag)


def deliver_rna(
    ignore_alignment,
    raw_folder,
    expression_folder,
    var_folder,
    alignment_folder,
    reports_folder,
    pcgr_folder,
    patient,
    out_base_path,
    in_base_path,
    in_base_path_abacus_rawdata,
    location_endpoint,
    file_dict,
    file_dict_rawdata,
    metrics_dict
    ):

    metric_mapping = {
        "raw_mean_coverage": "Raw_Mean_Coverage",
        "raw_duplication_rate": "Raw_Duplication_Rate",
        "raw_median_insert_size": "Raw_Median_Insert_Size",
        "raw_mean_insert_size": "Raw_Mean_Insert_Size",
        "raw_reads_count": "Raw_Reads_Count",
        "median_insert_size": "Median_Insert_Size",
        "mean_insert_size": "Mean_Insert_Size",
        "expression_profiling_efficiency": "WTS_Expression_Profiling_Efficiency",
        "aligned_reads_ratio": "WTS_Aligned_Reads",
        "rrna_rate": "WTS_rRNA_contamination",
    }

    # Compile the regex patterns
    expression_pattern = re.compile(r"abundance_(transcripts|genes)\.tsv$")
    variant_pattern = re.compile(r"\.hc\.vt\.annot(\.flt)?\.vcf\.gz")
    reports_pattern = re.compile(r"(\.multiqc|\.pcgr_acmg\.grch38\.flexdb)\.html$|\.putative_driver_fusions\.tsv$")
    pcgr_pattern = re.compile(r"(\.pcgr_acmg\.grch38\.(maf|snvs_indels\.tiers\.tsv))$")

    for sample in patient["sample"]:
        sample_name = sample["name"]
        for readset in sample["readset"]:
            readset_name = readset["name"]
            for file in readset["file"]:
                file_name = file["name"]
                # Handle Abacus rawdata path located in /lb/robot/research/freezeman-processing/novaseqx/
                if location_endpoint == "abacus" and file["location"].startswith(in_base_path_abacus_rawdata):
                    file_location = remove_path_parts(file["location"], in_base_path_abacus_rawdata)
                    # To workaround issue with RP naming regarding raw data being non unique we have to use file_location for file_name
                    file_dict_rawdata[file_location] = os.path.join(remove_path_parts(raw_folder, out_base_path), file_name)
                    continue
                # Usual case
                file_location = remove_path_parts(file["location"], in_base_path)
                # raw_data
                if "MAIN/raw_reads/" in file_location or "GQ_STAGING" in file_location:
                    file_dict[file_location] = os.path.join(remove_path_parts(raw_folder, out_base_path), file_name)
                # expression
                elif expression_pattern.search(file_name):
                    # Rename to include sample name in file
                    file_name = f"{sample_name}.{file_name}"
                    file_dict[file_location] = os.path.join(remove_path_parts(expression_folder, out_base_path), file_name)
                # variants
                elif variant_pattern.search(file_name):
                    # Rename as it's delivered as filt and not flt
                    file_name = file_name.replace("flt", "filt")
                    file_dict[file_location] = os.path.join(remove_path_parts(var_folder, out_base_path), file_name)
                # reports
                elif reports_pattern.search(file_name):
                    if "putative_driver_fusions" in file_name:
                        file_name = file_name.replace("putative_driver_fusions", "anno_fuse")
                    # Anno Fuse is not renamed by patient but keeping sample name
                    else:
                        # Rename to patient name instead of sample name
                        file_name = file_name.replace(sample_name, patient["name"])
                        match = reports_pattern.search(file_name)
                        if match:
                            # Insert '_R' before the matched part to differentiate between DNA and RNA reports
                            start = match.start()
                            file_name = file_name[:start] + "_R" + file_name[start:]
                        # Rename pcgr report
                        if "pcgr_acmg" in file_name:
                            file_name = file_name.replace("_acmg.grch38.flexdb", "")
                    file_dict[file_location] = os.path.join(remove_path_parts(reports_folder, out_base_path), file_name)
                # reports/pcgr
                elif pcgr_pattern.search(file_name):
                    match = pcgr_pattern.search(file_name)
                    if match:
                        # Insert '_D' before the matched part to differentiate between DNA and RNA reports
                        start = match.start()
                        file_name = file_name[:start] + "_R" + file_name[start:]
                    # Rename to patient name instead of sample name
                    file_name = file_name.replace(sample_name, patient["name"])
                    file_dict[file_location] = os.path.join(remove_path_parts(pcgr_folder, out_base_path), file_name)
                # alignment has to go as last if because otehr regex files are located under alignment folder
                elif "MAIN/alignment/" in file_location and not ignore_alignment:
                    # Rename to remove extra part from rna bams
                    file_name = file_name.replace("sorted.mdup.split.recal", "variants")
                    file_dict[file_location] = os.path.join(remove_path_parts(alignment_folder, out_base_path), file_name)
            for metric in readset["metric"]:
                metric_name = metric["name"]
                metric_value = metric["value"]
                metric_flag = metric["flag"]
                # metric_aggregate = metric["aggregate"]
                # update_metrics(metrics_dict, metric_mapping, sample_name, readset_name, metric_name, metric_value, metric_flag)
                update_metrics_dict(metrics_dict, metric_mapping, sample_name, readset_name, metric_name, metric_value, metric_flag)


def initialize_metrics(metric_mapping):
    metrics = {key: "NA" for key in metric_mapping.values()}
    metrics["Flags"] = "NA"
    metrics["Fails"] = "NA"
    return metrics


def update_metrics(existing_metrics, new_metrics_dict):
    updated_metrics = {row['Readset']: row for row in existing_metrics}
    for readset_name, metrics in new_metrics_dict.items():
        if readset_name in updated_metrics:
            updated_metrics[readset_name].update(metrics)
        else:
            new_row = {
                'Sample': metrics.get('sample_name', 'NA'),
                'Readset': readset_name,
                'WGS_Bases_Over_Q30': metrics.get('WGS_Bases_Over_Q30', 'NA'),
                'WGS_Min_Aligned_Reads_Delivered': metrics.get('WGS_Min_Aligned_Reads_Delivered', 'NA'),
                'Raw_Mean_Coverage': metrics.get('Raw_Mean_Coverage', 'NA'),
                'WGS_Dedup_Coverage': metrics.get('WGS_Dedup_Coverage', 'NA'),
                'Raw_Duplication_Rate': metrics.get('Raw_Duplication_Rate', 'NA'),
                'Raw_Median_Insert_Size': metrics.get('Raw_Median_Insert_Size', 'NA'),
                'Raw_Mean_Insert_Size': metrics.get('Raw_Mean_Insert_Size', 'NA'),
                'Median_Insert_Size': metrics.get('Median_Insert_Size', 'NA'),
                'Mean_Insert_Size': metrics.get('Mean_Insert_Size', 'NA'),
                'WGS_Contamination': metrics.get('WGS_Contamination', 'NA'),
                'Concordance': metrics.get('Concordance', 'NA'),
                'Purity': metrics.get('Purity', 'NA'),
                'Raw_Reads_Count': metrics.get('Raw_Reads_Count', 'NA'),
                'WTS_Exonic_Rate': metrics.get('WTS_Exonic_Rate', 'NA'),
                'WTS_Aligned_Reads': metrics.get('WTS_Aligned_Reads', 'NA'),
                'WTS_rRNA_contamination': metrics.get('WTS_rRNA_contamination', 'NA'),
                'WTS_Expression_Profiling_Efficiency': metrics.get('WTS_Expression_Profiling_Efficiency', 'NA'),
                'Flags': metrics.get('Flags', 'NA'),
                'Fails': metrics.get('Fails', 'NA')
            }
            updated_metrics[readset_name] = new_row
    return list(updated_metrics.values())

def update_metrics_dict(metrics_dict, metric_mapping, sample_name, readset_name, metric_name, metric_value, metric_flag):
    if readset_name not in metrics_dict:
        metrics_dict[readset_name] = initialize_metrics(metric_mapping)
        metrics_dict[readset_name]["Sample"] = sample_name

    if metric_name in metric_mapping:
        metrics_dict[readset_name][metric_mapping[metric_name]] = metric_value

        if metric_flag == "WARNING":
            # print(f"Sample name: {sample_name}, Metric name: {metric_name}, Metric flag: {metric_flag}")
            if metrics_dict[readset_name]["Flags"] == "NA":
                metrics_dict[readset_name]["Flags"] = metric_mapping[metric_name]
            else:
                metrics_dict[readset_name]["Flags"] += f";{metric_mapping[metric_name]}"
        elif metric_flag == "FAILED":
            if metrics_dict[readset_name]["Fails"] == "NA":
                metrics_dict[readset_name]["Fails"] = metric_mapping[metric_name]
            else:
                metrics_dict[readset_name]["Fails"] += f";{metric_mapping[metric_name]}"

def generate_key_metrics(metrics_dict):
    metrics_csv = []
    headers = "Sample,Readset,WGS_Bases_Over_Q30,WGS_Min_Aligned_Reads_Delivered,Raw_Mean_Coverage,WGS_Dedup_Coverage,Raw_Duplication_Rate,Raw_Median_Insert_Size,Raw_Mean_Insert_Size,Median_Insert_Size,Mean_Insert_Size,WGS_Contamination,Concordance,Purity,Raw_Reads_Count,WTS_Exonic_Rate,WTS_Aligned_Reads,WTS_rRNA_contamination,WTS_Expression_Profiling_Efficiency,Flags,Fails"
    metrics_csv.append(headers)
    for readset_name, metrics in metrics_dict.items():
        metrics_csv.append(",".join([
            metrics.get('sample_name', 'NA'),
            readset_name,
            metrics.get('WGS_Bases_Over_Q30', 'NA'),
            metrics.get('WGS_Min_Aligned_Reads_Delivered', 'NA'),
            metrics.get('Raw_Mean_Coverage', 'NA'),
            metrics.get('WGS_Dedup_Coverage', 'NA'),
            metrics.get('Raw_Duplication_Rate', 'NA'),
            metrics.get('Raw_Median_Insert_Size', 'NA'),
            metrics.get('Raw_Mean_Insert_Size', 'NA'),
            metrics.get('Median_Insert_Size', 'NA'),
            metrics.get('Mean_Insert_Size', 'NA'),
            metrics.get('WGS_Contamination', 'NA'),
            metrics.get('Concordance', 'NA'),
            metrics.get('Purity', 'NA'),
            metrics.get('Raw_Reads_Count', 'NA'),
            metrics.get('WTS_Exonic_Rate', 'NA'),
            metrics.get('WTS_Aligned_Reads', 'NA'),
            metrics.get('WTS_rRNA_contamination', 'NA'),
            metrics.get('WTS_Expression_Profiling_Efficiency', 'NA'),
            metrics.get('Flags', 'NA'),
            metrics.get('Fails', 'NA')
        ]))
    return "\n".join(metrics_csv)


def parse_existing_warnings(warnings_content):
    # If the content is empty or just whitespace, use a minimal fallback
    if not warnings_content.strip():
        warnings_content = "<table></table>"

    # Wrap the content in a full HTML structure to ensure soup.body exists
    wrapped_content = f"<html><head></head><body>{warnings_content}</body></html>"
    soup = BeautifulSoup(wrapped_content, 'html.parser')

    table = soup.find('table')
    warnings_dict = {}
    if table:
        rows = table.find_all('tr')[1:]  # Skip header row
        for row in rows:
            cols = row.find_all('td')
            if len(cols) >= 3:
                sample_name = cols[0].text.strip()
                flags = cols[1].text.strip()
                fails = cols[2].text.strip()
                for col in cols:
                    col.string = col.text.replace('\xa0', '').strip()
                warnings_dict[sample_name] = {
                    'Flags': flags,
                    'Fails': fails,
                    'row': row
                }
    return warnings_dict, soup


def update_warnings(existing_warnings_dict, metrics_dict, soup):
    for _, metrics in metrics_dict.items():
        sample_name = metrics.get('sample_name', 'NA')
        if sample_name in existing_warnings_dict:
            existing_warnings_dict[sample_name]['Flags'] = metrics.get('Flags', 'NA')
            existing_warnings_dict[sample_name]['Fails'] = metrics.get('Fails', 'NA')
            existing_warnings_dict[sample_name]['row'].find_all('td')[1].string = metrics.get('Flags', 'NA')
            existing_warnings_dict[sample_name]['row'].find_all('td')[2].string = metrics.get('Fails', 'NA')
        else:
            new_row = soup.new_tag('tr')
            new_row.append(soup.new_tag('td', style="text-align: left; padding: 8px; border: 1px solid #ddd;"))
            new_row.append(soup.new_tag('td', style="text-align: left; padding: 8px; border: 1px solid #ddd;"))
            new_row.append(soup.new_tag('td', style="text-align: left; padding: 8px; border: 1px solid #ddd;"))
            new_row.contents[0].string = sample_name
            new_row.contents[1].string = metrics.get('Flags', 'NA')
            new_row.contents[2].string = metrics.get('Fails', 'NA')
            existing_warnings_dict[sample_name] = {
                'Flags': metrics.get('Flags', 'NA'),
                'Fails': metrics.get('Fails', 'NA'),
                'row': new_row
            }
    return existing_warnings_dict


def generate_updated_warnings_content(warnings_dict, soup):
    table = soup.find('table')

    if not table:
        # Create a new table if it doesn't exist
        table = soup.new_tag('table')
        thead = soup.new_tag('thead')
        header_row = soup.new_tag('tr')
        for header in ['Sample Name', 'Flags', 'Fails']:
            th = soup.new_tag('th')
            th.string = header
            header_row.append(th)
        thead.append(header_row)
        table.append(thead)

        tbody = soup.new_tag('tbody')
        table.append(tbody)
        soup.body.append(table)  # Or insert it appropriately
    else:
        # Ensure thead exists
        thead = table.find('thead')
        if not thead:
            thead = soup.new_tag('thead')
            header_row = soup.new_tag('tr')
            for header in ['Sample Name', 'Flags', 'Fails']:
                th = soup.new_tag('th')
                th.string = header
                header_row.append(th)
            thead.append(header_row)
            table.insert(0, thead)  # Insert before tbody
        # Ensure tbody exists
        tbody = table.find('tbody')
        if not tbody:
            tbody = soup.new_tag('tbody')
            table.append(tbody)
        else:
            tbody.clear()

    for metrics in warnings_dict.values():
        tbody.append(metrics['row'])

    # Ensure the <head> element exists
    if not soup.head:
        head = soup.new_tag('head')
        soup.insert(0, head)
    else:
        head = soup.head

    # Add CSS styles for better readability
    style = soup.new_tag('style')
    style.string = """
    table {
        width: 100%;
        border-collapse: collapse;
        margin: 25px 0;
        font-size: 18px;
        text-align: left;
    }
    th, td {
        padding: 12px;
        border: 1px solid #ddd;
    }
    th {
        background-color: #f2f2f2;
    }
    tr:nth-child(even) {
        background-color: #f9f9f9;
    }
    tr:hover {
        background-color: #f1f1f1;
    }
    """
    head.append(style)

    return str(soup)



def file_exist_check(file, files_list):
    """Checking if the file exists in the list of files."""
    # Extract basename from the file path to check only file and not full path
    for i, (path, timestamp) in enumerate(files_list):
        if os.path.basename(path) == file:
            logger.debug(f"File {file} found with timestamp {timestamp}")
            # File found, removing from list to accelerate future searches
            files_list.pop(i)
            return f":white_check_mark: {timestamp}"

    return ":clock2:"


def generate_readme(
    patient,
    sample_name_dna_n,
    sample_name_dna_t,
    sample_name_rna,
    all_delivered_files):
    # Add timestamp
    timestamp = datetime.datetime.today().strftime("%Y/%m/%d")

    dna_n_raw = "\n    ".join(
        f"* `{os.path.basename(file[0])}` *Raw DNA reads for the Normal sample* {file_exist_check(os.path.basename(file[0]), all_delivered_files)}"
        for file in all_delivered_files if re.search(r"raw_data/.*DN.*\.bam", file[0])
    )
    if dna_n_raw:
        dna_n_raw = "\n    " + dna_n_raw
    dna_t_raw = "\n    ".join(
        f"* `{os.path.basename(file[0])}` *Raw DNA reads for the Tumor sample* {file_exist_check(os.path.basename(file[0]), all_delivered_files)}"
        for file in all_delivered_files if re.search(r"raw_data/.*DT.*\.bam", file[0])
    )
    if dna_t_raw:
        dna_t_raw = "\n    " + dna_t_raw
    rna_raw = "\n    ".join(
        f"* `{os.path.basename(file[0])}` *Raw RNA reads for the Tumor sample* {file_exist_check(os.path.basename(file[0]), all_delivered_files)}"
        for file in all_delivered_files if re.search(r"raw_data/.*RT.*\.fastq\.gz", file[0])
    )
    if rna_raw:
        rna_raw = "\n    " + rna_raw
    linx = "\n        ".join(
        f"* `{os.path.basename(file[0])}` {file_exist_check(os.path.basename(file[0]), all_delivered_files)}"
        for file in all_delivered_files if re.search(r"svariants/linx/.+", file[0])
    )
    if linx:
        linx = "\n        " + linx
    tumor_pair_inis = "\n    ".join(
        f"* `{os.path.basename(file[0])}` *Parameters used in the tumor pair analysis* {file_exist_check(os.path.basename(file[0]), all_delivered_files)}"
        for file in all_delivered_files if re.search(fr"parameters/{patient}\.TumourPair\..*ini", file[0])
    )
    if tumor_pair_inis:
        tumor_pair_inis = "\n    " + tumor_pair_inis
    rna_ligth_inis = "\n    ".join(
        f"* `{os.path.basename(file[0])}` *Parameters used in the RNA expression analysis* {file_exist_check(os.path.basename(file[0]), all_delivered_files)}"
        for file in all_delivered_files if re.search(fr"parameters/{patient}\.RNA\.Light\..*ini", file[0])
    )
    if rna_ligth_inis:
        rna_ligth_inis = "\n    " + rna_ligth_inis
    rna_cancer_inis = "\n    ".join(
        f"* `{os.path.basename(file[0])}` *Parameters used in the RNA variant analysis* {file_exist_check(os.path.basename(file[0]), all_delivered_files)}"
        for file in all_delivered_files if re.search(fr"parameters/{patient}\.RNA\.Variants\..*ini", file[0])
    )
    if rna_cancer_inis:
        rna_cancer_inis = "\n    " + rna_cancer_inis
    data = f"""This directory contains the delivered data for **{patient}** processed by the Canadian Centre for Computational Genomics.
The data will be updated as it becomes available and as such many files may be missing from RNA or DNA upon initial creation of this directory
Should you have concerns, questions, or suggestions, please contact the analysis team at moh-q@computationalgenomics.ca

Within this directory you will find the results of the analysis for a single patient; if multiple bam files are present in the `raw_data` folder it means TopTup(s) have been necessary.
When :white_check_mark: is present the file is available and the date at which it has been delivered is written next to it, when :clock2: is present the file is not ready for delivery yet.

* `Readme.html` *This file* :white_check_mark:
* [`Warnings.html`](Warnings.html) *Contains details of any warnings and whether they caused a failure of this analysis* {file_exist_check("Warnings.html", all_delivered_files)}
* [`Methods.html`](Methods.html) *Contains details on pipelines and references used for the analysis* {file_exist_check("Methods.html", all_delivered_files)}
* `Key_metrics.csv` *File with metrics for the patient in csv format* {file_exist_check("Key_metrics.csv", all_delivered_files)}
* `raw_data/` *Contains all of the BAMs/FASTQs from the sequencer. BAM files here include both mapped and unmapped reads and can be converted to the FASTQ format with tools such as SamToFastq.*{dna_n_raw}{dna_t_raw}{rna_raw}
* `variants/` *Contains the vcfs related to variant calls*
    * `{patient}.ensemble.germline.vt.annot.vcf.gz` *Germline Variants found in any of the callers for DNA Sample* {file_exist_check(f"{patient}.ensemble.germline.vt.annot.vcf.gz", all_delivered_files)}
    * `{patient}.ensemble.germline.vt.annot.vcf.gz.tbi` *Index of Germline Variants found in any of the callers for DNA Sample* {file_exist_check(f"{patient}.ensemble.germline.vt.annot.vcf.gz.tbi", all_delivered_files)}
    * `{patient}.ensemble.somatic.vt.annot.vcf.gz` *Somatic Variants found in any of the callers for DNA Sample* {file_exist_check(f"{patient}.ensemble.somatic.vt.annot.vcf.gz", all_delivered_files)}
    * `{patient}.ensemble.somatic.vt.annot.vcf.gz.tbi` *Index of Somatic Variants found in any of the callers for DNA Sample* {file_exist_check(f"{patient}.ensemble.somatic.vt.annot.vcf.gz.tbi", all_delivered_files)}
    * `{patient}.hc.vt.annot.vcf.gz` *Variants found using RNA sample* {file_exist_check(f"{patient}.hc.vt.annot.vcf.gz", all_delivered_files)}
    * `{patient}.hc.vt.annot.vcf.gz.tbi` *Index of Variants found using RNA sample* {file_exist_check(f"{patient}.hc.vt.annot.vcf.gz.tbi", all_delivered_files)}
    * `{patient}.hc.vt.annot.filt.vcf.gz` *Variants found using RNA sample; annotated with RNAEdits and filtered with coverage >=10x and VAF >=5%* {file_exist_check(f"{patient}.hc.vt.annot.filt.vcf.gz", all_delivered_files)}
    * `{patient}.hc.vt.annot.filt.vcf.gz.tbi` *Index of Variants found using RNA sample; annotated with RNAEdits and filtered with coverage >=10x and VAF >=5%* {file_exist_check(f"{patient}.hc.vt.annot.filt.vcf.gz.tbi", all_delivered_files)}
    * `caller_vcfs/` *Contains the vcfs produced from individual callers on the DNA samples*
        * `{patient}.mutect2.somatic.vt.vcf.gz` *Somatic results for mutect2* {file_exist_check(f"{patient}.mutect2.somatic.vt.vcf.gz", all_delivered_files)}
        * `{patient}.strelka2.germline.vt.vcf.gz` *Germline results for strelka2* {file_exist_check(f"{patient}.strelka2.germline.vt.vcf.gz", all_delivered_files)}
        * `{patient}.strelka2.somatic.vt.vcf.gz` *Somatic results for strelka2* {file_exist_check(f"{patient}.strelka2.somatic.vt.vcf.gz", all_delivered_files)}
        * `{patient}.vardict.germline.vt.vcf.gz` *Germline results for vardict* {file_exist_check(f"{patient}.vardict.germline.vt.vcf.gz", all_delivered_files)}
        * `{patient}.vardict.somatic.vt.vcf.gz` *Somatic results for vardict* {file_exist_check(f"{patient}.vardict.somatic.vt.vcf.gz", all_delivered_files)}
        * `{patient}.varscan2.germline.vt.vcf.gz` *Germline results for varscan2* {file_exist_check(f"{patient}.varscan2.germline.vt.vcf.gz", all_delivered_files)}
        * `{patient}.varscan2.somatic.vt.vcf.gz` *Somatic results for varscan2* {file_exist_check(f"{patient}.varscan2.somatic.vt.vcf.gz", all_delivered_files)}
* `svariants/` *Contains the files related to structural variants*
    * `{patient}.gridss.vcf.gz` *Unfiltered structural variant calls with GRIDSS* {file_exist_check(f"{patient}.gridss.vcf.gz", all_delivered_files)}
    * `{patient}.gripss.filtered.somatic.vcf.gz` *Annotated and filtered somatic structural variant calls from GRIDSS using GRIPSS* {file_exist_check(f"{patient}.gripss.filtered.somatic.vcf.gz", all_delivered_files)}
    * `{patient}.gripss.filtered.germline.vcf.gz` *Annotated and filtered germline structural variant calls from GRIDSS using GRIPSS* {file_exist_check(f"{patient}.gripss.filtered.germline.vcf.gz", all_delivered_files)}
    * `{patient}.driver.catalog.somatic.tsv` *Driver somatic structural variant calls* {file_exist_check(f"{patient}.driver.catalog.somatic.tsv", all_delivered_files)}
    * `{patient}.driver.catalog.germline.tsv` *Driver germline structural variant calls* {file_exist_check(f"{patient}.driver.catalog.germline.tsv", all_delivered_files)}
    * `{patient}.circos.png` *Circos plot of all variants. Cf. https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md#circos* {file_exist_check(f"{patient}.circos.png", all_delivered_files)}
    * `{patient}.purple_ensemble.zip` {file_exist_check(f"{patient}.purple_ensemble.zip", all_delivered_files)}
    * `{patient}.purple_sv.zip` {file_exist_check(f"{patient}.purple_sv.zip", all_delivered_files)}
    * `linx/` *Contains structural variant annotated and visualized by LINX Cf. https://github.com/hartwigmedical/hmftools/tree/master/linx*{linx}
* `raw_cnv/` *Contains the raw copy number calls for each patient DNA*
    * `{patient}.cnvkit.vcf.gz` *Raw cnvkit output* {file_exist_check(f"{patient}.cnvkit.vcf.gz", all_delivered_files)}
    * `{patient}.cnvkit.vcf.gz.tbi` *Index of Raw cnvkit output* {file_exist_check(f"{patient}.cnvkit.vcf.gz.tbi", all_delivered_files)}
* `alignment/` *Contains the alignment data for each sample*
    * `{sample_name_dna_n}.bam` *Alignment of normal against the reference* {file_exist_check(f"{sample_name_dna_n}.bam", all_delivered_files)}
    * `{sample_name_dna_n}.bam.bai` *Index of Alignment of normal against the reference* {file_exist_check(f"{sample_name_dna_n}.bam.bai", all_delivered_files)}
    * `{sample_name_dna_t}.bam` *Alignment of tumor against the reference* {file_exist_check(f"{sample_name_dna_t}.bam", all_delivered_files)}
    * `{sample_name_dna_t}.bam.bai` *Index of Alignment of tumor against the reference* {file_exist_check(f"{sample_name_dna_t}.bam.bai", all_delivered_files)}
    * `{sample_name_rna}.variants.bam` *Alignment of tumor RNA against the reference used in variants analysis* {file_exist_check(f"{sample_name_rna}.variants.bam", all_delivered_files)}
    * `{sample_name_rna}.variants.bam.bai` *Index of Alignment of tumor RNA against the reference used in variants analysis* {file_exist_check(f"{sample_name_rna}.variants.bam.bai", all_delivered_files)}
* `expression/` *Contains the transcripts and genes abundance estimation from Kallisto*
    * `{sample_name_rna}.abundance_transcripts.tsv` *Table with transcript abundance from Kallisto* {file_exist_check(f"{sample_name_rna}.abundance_transcripts.tsv", all_delivered_files)}
    * `{sample_name_rna}.abundance_genes.tsv` *Table with gene abundance from Kallisto* {file_exist_check(f"{sample_name_rna}.abundance_genes.tsv", all_delivered_files)}
* `reports/` *Contains the reports for the experiment*
    * [`{patient}_D.multiqc.html`](reports/{patient}_D.multiqc.html) *QC report for the DNA analysis* {file_exist_check(f"{patient}_D.multiqc.html", all_delivered_files)}
    * [`{patient}_R.multiqc.html`](reports/{patient}_R.multiqc.html) *QC report for the RNA analysis* {file_exist_check(f"{patient}_R.multiqc.html", all_delivered_files)}
    * `{sample_name_rna}.anno_fuse.tsv` *TSV for fusions detected using RN, to be loaded into Excel sheet for easier display* {file_exist_check(f"{sample_name_rna}.anno_fuse.tsv", all_delivered_files)}
    * [`{patient}_D.pcgr.html`](reports/{patient}_D.pcgr.html) *Personal Cancer Genome Reporter report for the DNA analysis; Cf. https://pcgr.readthedocs.io/en/latest* {file_exist_check(f"{patient}_D.pcgr.html", all_delivered_files)}
    * [`{patient}_R.pcgr.html`](reports/{patient}_R.pcgr.html) *Personal Cancer Genome Reporter report for the RNA analysis; Cf. https://pcgr.readthedocs.io/en/latest* {file_exist_check(f"{patient}_R.pcgr.html", all_delivered_files)}
    * `{patient}_D.cpsr.zip` *Cancer Predisposition Sequencing Reporter archive for the DNA analysis; Cf. https://sigven.github.io/cpsr/* {file_exist_check(f"{patient}_D.cpsr.zip", all_delivered_files)}
    * `pcgr/` *Contains raw tables used to generate PCGR reports*
        * `{patient}_D.acmg.grch38.maf` {file_exist_check(f"{patient}_D.acmg.grch38.maf", all_delivered_files)}
        * `{patient}_D.acmg.grch38.snvs_indels.tiers.tsv` {file_exist_check(f"{patient}_D.acmg.grch38.snvs_indels.tiers.tsv", all_delivered_files)}
        * `{patient}_D.acmg.grch38.cna_segments.tsv.gz` {file_exist_check(f"{patient}_D.acmg.grch38.cna_segments.tsv.gz", all_delivered_files)}
        * `{patient}_R.acmg.grch38.maf` {file_exist_check(f"{patient}_R.acmg.grch38.maf", all_delivered_files)}
        * `{patient}_R.acmg.grch38.snvs_indels.tiers.tsv` {file_exist_check(f"{patient}_R.acmg.grch38.snvs_indels.tiers.tsv", all_delivered_files)}
* `parameters/` *Contains the records of all the Parameters used in the pipeline analysis*{tumor_pair_inis}{rna_ligth_inis}{rna_cancer_inis}

Generated {timestamp}."""
    html = markdown.markdown(data, extensions=extensions, extension_configs=extension_configs)
    return html

def generate_methods():
    # get_local_file_log(methods_file, log, updated, old_log)
    data = r"""# Methods
Using an Illumina NovaSeq 6000 instrument, Whole Genome Sequencing (WGS) was performed on the tumor and matched normal samples, with a target depth of coverage of 80X and 30X respectively. Similarly, Whole Transcriptome Sequencing (WTS) was done on tumour samples, with a target of 100 million paired-end reads per sample.

Bioinformatics analyses were performed using the [GenPipes][GenPipes_BB][^GenPipes_] Tumor-Pair and RNA-seq analytical pipelines (detailed documentation can be found [here][GenPipes_RTD]). Specific parameter values and reference databases are tracked in corresponding *\*.ini* files found under the `parameters` directory. An explanation of the various steps is detailed below.

## WGS
For WGS samples, [GATK][GATK_]’s best practices and procedures were followed. Quality trimmed and adapter-clipped reads were first aligned to the GRCh38 reference[^genome_ref_] with BWA-MEM[^BWA-MEM_]. Alignments were then sorted, realigned around Indels, and marked for duplicates. Base qualities were improved using Base Quality Score Recalibration (BQSR). A MultiQC[^MultiQC_] report was generated per patient to flag any inconsistency in overall coverage, QC bias, tumor purity ([PURPLE][purple_]) and normal or tumor contamination and concordance estimations (ConPair[^ConPair_]).  Somatic and germline calls were generated using an ensemble approach combining four independent variant callers: GATK MuTect2, Strelka2[^Strelka2_], VarDict[^VarDict_] and VarScan2[^VarScan2_]. Somatic and germline variants identified in two or more callers were further annotated and prioritized using the PCGR/CSPR[^PCGR_] reporting system to classify somatic and germline calls using ACMP/AMP classification, perform tumor mutation burden (TMB) estimation, microsatellite instability (MSI) classification, mutational signature estimations, and kataegis detection. Structural Variants were called using GRIDSS with PURPLE and LINX[^GRIDSS_PURPLE_LINX_] as an interpretation tool, and CNAs were called with CNVKit[^CNVKit_]. 

## WTS
For WTS, transcript abundance was estimated using the Kallisto[^Kallisto_] pseudoaligner from quality trimmed and adapter-clipped reads. For WTS variant calling, a full alignment to the same GRCh38 reference was performed using STAR[^STAR_] and the alignments were sorted, realigned and duplicates were marked using GATK best practices and procedures. Fusions in WTS data were assessed using both STAR-FUSION[^STAR-FUSION_] and Arriba[^Arriba_] fusion callers and reported using AnnoFuse[^AnnoFuse_]. Finally, aligned reads were used to call variants with the GATK haplotype caller. Variant calls were then filtered, annotated and prioritized using PCGR/CSPR reporting systems.

In both WGS and WTS, reports were generated with MultiQC for visualization of key quality metrics and allow for manual validation.    

# References
[^GenPipes_]: Bourgey M, Dali R, Eveleigh R, Chen KC, Letourneau L, Fillon J, Michaud M, Caron M, Sandoval J, Lefebvre F, Leveque G, Mercier E, Bujold D, Marquis P, Van PT, Anderson de Lima Morais D, Tremblay J, Shao X, Henrion E, Gonzalez E, Quirion PO, Caron B, Bourque G. GenPipes: an open-source framework for distributed and scalable genomic analyses. Gigascience. 2019 Jun 1;8(6):giz037. doi: 10.1093/gigascience/giz037. PMID: 31185495; PMCID: PMC6559338.

[GenPipes_RTD]: https://genpipes.readthedocs.io


[GenPipes_BB]: https://bitbucket.org/mugqic/genpipes

[GATK_]: https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows

[^genome_ref_]: GCA_000001405.15 no alt analysis set downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz


[^BWA-MEM_]: Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN]


[^MultiQC_]: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PMID: 27312411; PMCID: PMC5039924.


[purple_]: https://github.com/hartwigmedical/hmftools/tree/master/purple

[^ConPair_]: Bergmann EA, Chen BJ, Arora K, Vacic V, Zody MC. Conpair: concordance and contamination estimator for matched tumor-normal pairs. Bioinformatics. 2016 Oct 15;32(20):3196-3198. doi: 10.1093/bioinformatics/btw389. Epub 2016 Jun 26. PMID: 27354699; PMCID: PMC5048070.

[^Strelka2_]: Kim S, Scheffler K, Halpern AL, Bekritsky MA, Noh E, Källberg M, Chen X, Kim Y, Beyter D, Krusche P, Saunders CT. Strelka2: fast and accurate calling of germline and somatic variants. Nat Methods. 2018 Aug;15(8):591-594. doi: 10.1038/s41592-018-0051-x. Epub 2018 Jul 16. PMID: 30013048.

[^VarDict_]: Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, Johnson J, Dougherty B, Barrett JC, Dry JR. VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. Nucleic Acids Res. 2016 Jun 20;44(11):e108. doi: 10.1093/nar/gkw227. Epub 2016 Apr 7. PMID: 27060149; PMCID: PMC4914105.

[^VarScan2_]: Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, Miller CA, Mardis ER, Ding L, Wilson RK. VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Res. 2012 Mar;22(3):568-76. doi: 10.1101/gr.129684.111. Epub 2012 Feb 2. PMID: 22300766; PMCID: PMC3290792.

[^PCGR_]: Nakken S, Fournous G, Vodák D, Aasheim LB, Myklebost O, Hovig E. Personal Cancer Genome Reporter: variant interpretation report for precision oncology. Bioinformatics. 2018 May 15;34(10):1778-1780. doi: 10.1093/bioinformatics/btx817. PMID: 29272339; PMCID: PMC5946881.

[^GRIDSS_PURPLE_LINX_]: GRIDSS, PURPLE, LINX: Unscrambling the tumor genome via integrated analysis of structural variation and copy number. Daniel L. Cameron, Jonathan Baber, Charles Shale, Anthony T. Papenfuss, Jose Espejo Valle-Inclan, Nicolle Besselink, Edwin Cuppen, Peter Priestley. bioRxiv 781013; doi: https://doi.org/10.1101/781013

[^CNVKit_]: Talevich E, Shain AH, Botton T, Bastian BC. CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. PLoS Comput Biol. 2016 Apr 21;12(4):e1004873. doi: 10.1371/journal.pcbi.1004873. PMID: 27100738; PMCID: PMC4839673.

[^Kallisto_]: Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol. 2016 May;34(5):525-7. doi: 10.1038/nbt.3519. Epub 2016 Apr 4. Erratum in: Nat Biotechnol. 2016 Aug 9;34(8):888. PMID: 27043002.

[^STAR_]: Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635. Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.

[^STAR-FUSION_]: STAR-Fusion: Fast and Accurate Fusion Transcript Detection from RNA-Seq. Brian J. Haas, Alex Dobin, Nicolas Stransky, Bo Li, Xiao Yang, Timothy Tickle, Asma Bankapur, Carrie Ganote, Thomas G. Doak, Nathalie Pochet, Jing Sun, Catherine J. Wu, Thomas R. Gingeras, Aviv Regev. bioRxiv 120295; doi: https://doi.org/10.1101/120295

[^Arriba_]: Uhrig S, Ellermann J, Walther T, Burkhardt P, Fröhlich M, Hutter B, Toprak UH, Neumann O, Stenzinger A, Scholl C, Fröhling S, Brors B. Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Res. 2021 Mar;31(3):448-460. doi: 10.1101/gr.257246.119. Epub 2021 Jan 13. PMID: 33441414; PMCID: PMC7919457.

[^AnnoFuse_]: Gaonkar KS, Marini F, Rathi KS, Jain P, Zhu Y, Chimicles NA, Brown MA, Naqvi AS, Zhang B, Storm PB, Maris JM, Raman P, Resnick AC, Strauch K, Taroni JN, Rokita JL. annoFuse: an R Package to annotate, prioritize, and interactively explore putative oncogenic RNA fusions. BMC Bioinformatics. 2020 Dec 14;21(1):577. doi: 10.1186/s12859-020-03922-7. PMID: 33317447; PMCID: PMC7737294."""
    html = markdown.markdown(data, extensions=extensions, extension_configs=extension_configs)
    # with open(methods_file, "w", encoding="utf-8", errors="xmlcharrefreplace") as out_file:
    #     out_file.write(html)
    return html

def read_env_file(file_path):
    """Function to read environment variables from a file"""
    with open(file_path) as f:
        for line in f:
            if line.startswith('export'):
                key, value = line.split('=', 1)[0].split()[1], line.split('=', 1)[1].strip().strip('"')
                os.environ[key] = value

@retry(globus_sdk.GlobusAPIError, tries=10, delay=2, backoff=2)
def fetch_file_data(endpoint_id, path, transfer_client):
    """Function to fetch file data with retry mechanism"""
    response = transfer_client.operation_ls(endpoint_id, path=path)
    return response

def get_remote_file_mod_date(endpoint_id, file_path, transfer_client):
    """Function to get file modification date"""
    response = fetch_file_data(endpoint_id, file_path, transfer_client)
    for item in response['DATA']:
        if item['type'] == 'file' and item['name'] == file_path.split('/')[-1]:
            return parse(item['last_modified'])
    return None  # Return None if the file is not found

def transfer_files_with_sync(endpoint_id_src, endpoint_id_dest, file_dict, transfer_client, transfer_label, sync_level=2, encrypt=True):
    """Function to transfer files with sync level and track transferred/skipped files"""
    transfer_data = globus_sdk.TransferData(
        transfer_client,
        endpoint_id_src,
        endpoint_id_dest,
        label=transfer_label,
        sync_level=sync_level,
        encrypt_data=encrypt
        # recursive_symlinks="copy"
    )

    for src_file, dest_file in file_dict.items():
        transfer_data.add_item(src_file, dest_file)

    transfer_task = transfer_client.submit_transfer(transfer_data)
    print(f"Track your transfer here: https://app.globus.org/activity/{transfer_task['task_id']}")
    return transfer_task['task_id']

def get_transfer_event_log(transfer_client, task_id, s3_client, bucket_name):
    """Function to get the event log of a transfer task with timestamps"""
    events = transfer_client.task_event_list(task_id)
    transferred_files = []
    # skipped_files = []
    fault_events = []

    # Get fault events as they occur
    events = transfer_client.task_event_list(task_id)
    for event in events:
        if event.get('is_error', False):
            try:
                details = json.loads(event['details'])  # Parse the JSON string
                error_body = details.get('error', {}).get('body')
                if error_body:
                    fault_events.append(error_body)
                else:
                    fault_events.append(event.get('description', 'No description available'))
            except json.JSONDecodeError:
                fault_events.append(event.get('description', 'No description available'))


    # Check task status
    task = transfer_client.get_task(task_id)
    if task['status'] != 'SUCCEEDED':
        return transferred_files, fault_events

    # Get successful transfers
    success_transfers = transfer_client.task_successful_transfers(task_id)
    # print(success_transfers)
    for transfer in success_transfers:
        destination_path = transfer['destination_path'].replace("/~/", "")
        timestamp = transfer.get('timestamp')
        if not timestamp:
            # Fetch the timestamp from S3
            response = s3_client.head_object(Bucket=bucket_name, Key=destination_path)
            timestamp = response['LastModified'].isoformat()

        transferred_files.append((destination_path, timestamp))

    # Get skipped files
    # skipped_transfers = transfer_client.task_skipped_errors(task_id)
    # print(skipped_transfers)
    # for skip in skipped_transfers:
    #     skipped_files.append((skip['destination_path'], skip['timestamp']))

    return transferred_files, fault_events


def display_transfer_status(transfer_client, task_id, s3_client, bucket_name):
    """Function to display live transfer status and speed"""

    def signal_handler(sig, frame):
        logger.debug("Transfer interrupted. Fetching final event log...")
        stop_transfer(transfer_client, task_id)
        transferred_files, fault_events = get_transfer_event_log(transfer_client, task_id, s3_client, bucket_name)
        if transferred_files:
            logger.debug(f"Transferred files:\n {'\n '.join([file[0] for file in transferred_files])}")
        # if skipped_files:
        #     logger.info(f"Skipped files: {skipped_files}")
        if fault_events:
            logger.info(f"Fault events:\n {'\n '.join(list(fault_events))}")
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)

    # Wait until start_time is available
    start_time = None
    while not start_time:
        task = transfer_client.get_task(task_id)
        start_time = task.get('request_time', None)
        if not start_time:
            time.sleep(10)  # Wait for 10 seconds before checking again

    start_time = datetime.datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S%z")

    previous_message_length = 0

    while True:
        task = transfer_client.get_task(task_id)
        status = task['status']
        bytes_transferred = task.get('bytes_transferred', None)
        transfer_rate = task.get('effective_bytes_per_second', None)
        end_time = task.get('completion_time', None)

        # logger.info(f"Status: {status}")
        # logger.info(f"Transferred: {human_readable_size(bytes_transferred)}")
        # logger.info(f"Transfer rate: {human_readable_size(transfer_rate)}/second")

        # Calculate the elapsed time
        current_time = datetime.datetime.now(datetime.timezone.utc)
        if end_time:
            end_time = datetime.datetime.strptime(end_time, "%Y-%m-%dT%H:%M:%S%z")
            elapsed_time = end_time - start_time
        else:
            elapsed_time = current_time - start_time

            # Ensure elapsed_time is not negative
            if elapsed_time.total_seconds() < 0:
                elapsed_time = datetime.timedelta(seconds=0)

        duration = str(elapsed_time).split('.', maxsplit=1)[0] # Remove microseconds

        # Create the status message
        status_message = (
            f"Status: {status} | "
            f"Transferred: {human_readable_size(bytes_transferred)} | "
            f"Transfer rate: {human_readable_size(transfer_rate)}/second | "
            f"Duration: {duration}"
        )

        # Pad the message with spaces to clear any remaining characters
        padded_message = status_message.ljust(previous_message_length)
        previous_message_length = len(status_message)


        # Print the status message without a new line and flush the output
        sys.stdout.write(f"\r{padded_message}")
        sys.stdout.flush()


        transferred_files, fault_events = get_transfer_event_log(transfer_client, task_id, s3_client, bucket_name)
        # if skipped_files:
        #     logger.info(f"Skipped files: {skipped_files}")
        if fault_events:
            # Create the status message
            status_message = (
                f"Status: {status} | "
                f"Transferred: {human_readable_size(bytes_transferred)} | "
                f"Transfer rate: {human_readable_size(transfer_rate)}/second | "
                f"Duration: {duration} | "
                f"Fault events:\n {'\n '.join(list(fault_events))}"
            )

            # Pad the message with spaces to clear any remaining characters
            padded_message = status_message.ljust(previous_message_length)
            previous_message_length = len(status_message)


            # Print the status message without a new line and flush the output
            sys.stdout.write(f"\r{padded_message}")
            sys.stdout.flush()
            # logger.info(f"Fault events: {fault_events}")

        if status in ['SUCCEEDED', 'FAILED', 'CANCELED']:
            # To have a new line after the sys.stdout write and flush
            print()
            if transferred_files:
                logger.debug(f"Transferred files:\n {'\n '.join([file[0] for file in transferred_files])}")
            break

        time.sleep(60)


def stop_transfer(transfer_client, task_id):
    """Function to stop the transfer"""
    transfer_client.cancel_task(task_id)
    logger.info(f"\nTransfer task {task_id} has been canceled.")


def human_readable_size(size, decimal_places=2):
    """Convert bytes to a human-readable format."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024:
            return f"{size:.{decimal_places}f} {unit}"
        size /= 1024


def format_timestamps(file_list):
    """Function to format timestamps in a list of files"""
    formatted_files = [(file, datetime.datetime.strptime(timestamp, "%Y-%m-%dT%H:%M:%S%z").strftime("%Y/%m/%d")) for file, timestamp in file_list]
    return formatted_files


def parse_csv(file_content):
    csv_reader = csv.DictReader(StringIO(file_content))
    return list(csv_reader)

def generate_csv_content(metrics_list):
    headers = "Sample,Readset,WGS_Bases_Over_Q30,WGS_Min_Aligned_Reads_Delivered,Raw_Mean_Coverage,WGS_Dedup_Coverage,Raw_Duplication_Rate,Raw_Median_Insert_Size,Raw_Mean_Insert_Size,Median_Insert_Size,Mean_Insert_Size,WGS_Contamination,Concordance,Purity,Raw_Reads_Count,WTS_Exonic_Rate,WTS_Aligned_Reads,WTS_rRNA_contamination,WTS_Expression_Profiling_Efficiency,Flags,Fails"
    output = StringIO()
    csv_writer = csv.DictWriter(output, fieldnames=headers.split(','))
    csv_writer.writeheader()
    csv_writer.writerows(metrics_list)
    return output.getvalue()

def file_exists(s3_client, bucket_name, key):
    try:
        s3_client.head_object(Bucket=bucket_name, Key=key)
        return True
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == '404':
            return False
        raise e

def remove_path_parts(original_path, part_to_remove):
    """Function to remove a part of the path"""
    try:
        # Create Path object
        path_obj = Path(original_path)

        # Remove specified part of the path
        new_path = path_obj.relative_to(part_to_remove)

        # Convert Path object to string
        return str(new_path)
    except ValueError:
        # Handle case where part_to_remove is not in original_path
        return original_path

def list_s3_files(s3_client, bucket_name, prefix):
    """
    Function to list all files and their formatted timestamps in a given S3 bucket path
    Arguments:
        s3_client -- boto3 S3 client
        bucket_name -- name of the S3 bucket
        prefix -- prefix path in the S3 bucket
    Returns:
        files -- list of tuples containing file names and their formatted timestamps
    """
    files = []
    continuation_token = None
    retries = 5

    while True:
        try:
            if continuation_token:
                response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, ContinuationToken=continuation_token)
            else:
                response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix)

            files.extend([(item['Key'], item['LastModified'].strftime("%Y/%m/%d")) for item in response.get('Contents', [])])

            if response.get('IsTruncated'):
                continuation_token = response.get('NextContinuationToken')
            else:
                break
        except Exception as e:
            if retries > 0:
                retries -= 1
                time.sleep(2 ** (5 - retries))  # Exponential backoff
            else:
                raise e

    return files


if __name__ == '__main__':
    main()
