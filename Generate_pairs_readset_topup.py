#!/usr/bin/env python3

import glob
import sys
import os
import re
import datetime
import shutil
import argparse
import csv

from  DB_OPS import create_connection, Update_Samples_Table

def main():
    parser = argparse.ArgumentParser(prog='generate_genpipes_inputs.py', description="Generates readset file for GenPipes and move raw_reads in MAIN?raw_reads folder for MOH project processing.")
    parser.add_argument('--type', required=True, help="Type of analysis: either RNA or DNA", choices=['DNA', 'RNA'])
    parser.add_argument('--dry_run', required=False, help="Will not move raw_data reads generate GenPipes readset/pair file(s) (Default: False)", default=False, action='store_true')
    parser.add_argument('--force', required=False, help="If used the duplicates samples will be overwritten (Default: False)", default=False, action='store_true')
    parser.add_argument('--output_folder', required=False, help=f"Folder in which readsets will be written (Default: {os.getcwd()})", default=os.getcwd())
    args = parser.parse_args()

    sequencing_type = args.type

    beluga_moh_folder = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING"
    beluga_main_folder = os.path.join(beluga_moh_folder, "MAIN")
    beluga_db = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db"
    beluga_main_raw_reads_folder = os.path.join(beluga_main_folder, "raw_reads")
    beluga_transferred_raw_reads_folder = os.path.join(beluga_moh_folder, "raw_reads")

    date = datetime.datetime.today()
    date_formatted = date.strftime("%Y-%m-%dT%H.%M.%S")
    #Connect to the db
    connection = create_connection(beluga_db)

    readset_header = 'Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM'

    # Populate the lists
    rna_samples_list = []
    patient_dict = {}
    dna_patient_pair_dict = {}
    transferred_samples_list = os.listdir(beluga_transferred_raw_reads_folder)
    duplicates = []
    for sample in transferred_samples_list:
        if sample.startswith("MoHQ"):
            sample_object = Sample(sample)
            if sequencing_type == 'RNA' and sample.endswith('RT'):
                main_fastq = glob.glob(os.path.join(beluga_main_raw_reads_folder, sample, "*.fastq.gz"))
                main_fastq = [os.path.basename(file) for file in main_fastq]
                raw_fastq = glob.glob(os.path.join(beluga_transferred_raw_reads_folder, sample, "*.fastq.gz"))
                raw_fastq = [os.path.basename(file) for file in raw_fastq]
                if intersection(raw_fastq, main_fastq):
                    duplicates.append(sample)
                # if os.path.exists(os.path.join(beluga_main_raw_reads_folder, sample)):
                #     duplicates.append(sample)
                # rna_samples_list.append(sample)
                try:
                    patient_dict[sample_object.patient][sample_object.type] = sample_object
                    # dna_patient_pair_dict[sample_object.patient] = patient_dict[sample_object.patient]
                except KeyError:
                    patient_dict[sample_object.patient] = {sample_object.type: sample_object}
            elif sequencing_type == 'DNA' and (sample.endswith('DT') or sample.endswith('DN')):
                main_bam = glob.glob(os.path.join(beluga_main_raw_reads_folder, sample, "*.bam"))
                main_bam = [os.path.basename(file) for file in main_bam]
                raw_bam = glob.glob(os.path.join(beluga_transferred_raw_reads_folder, sample, "*.bam"))
                raw_bam = [os.path.basename(file) for file in raw_bam]
                if intersection(raw_bam, main_bam):
                    duplicates.append(sample)
                # if os.path.exists(os.path.join(beluga_main_raw_reads_folder, sample)):
                #     duplicates.append(sample)
                # sample_object = Sample(sample)
                try:
                    patient_dict[sample_object.patient][sample_object.type] = sample_object
                    dna_patient_pair_dict[sample_object.patient] = patient_dict[sample_object.patient]
                except KeyError:
                    patient_dict[sample_object.patient] = {sample_object.type: sample_object}

    if duplicates and not args.force:
        duplicates_file = os.path.join(args.output_folder, f"{date_formatted}_duplicates.txt")
        os.makedirs(os.path.dirname(duplicates_file), exist_ok=True)
        with open(duplicates_file, "w+", encoding="utf-8") as filename:
            for line in duplicates:
                filename.write(f"{line}\n")
        sys.exit(f"Samples are already present in in MAIN/raw_reads folder. Either use --force option to overwritte them or manually delete them. Cf. file {duplicates_file} for a list of duplicates.\nExiting...")

    # RNA
    if sequencing_type == 'RNA':
        if patient_dict:
            readset_rna_out = []
            # Check to see if any files are already present in the final directory
            # duplicate_samples(beluga_main_folder, beluga_main_raw_reads_folder, rna_samples_list, date_formatted)
            # readset generation
            for patient, sample in patient_dict.items():
                rna_sample = sample['RT'].sample
                cohort = sample['RT'].cohort
                institution = sample['RT'].institution
                transferred_readsets = glob.glob(os.path.join(beluga_transferred_raw_reads_folder, rna_sample, "*_readset.tsv"))
                if not transferred_readsets:
                    sys.exit(f"Sample {rna_sample} doesn't have a readset file, there is an issue with the transfer. Please check {os.path.join(beluga_transferred_raw_reads_folder, rna_sample)}\nExiting...")
                for transferred_readset in transferred_readsets:
                    with open(transferred_readset, 'rt') as readset_in:
                        reader = csv.reader(readset_in, delimiter="\t")
                        next(reader, None)
                        for line in reader:
                            # Making fastq/bam paths absolute
                            if line[-1]:
                                line[-1] = os.path.join(beluga_main_folder, line[-1])
                            if line[-2]:
                                line[-2] = os.path.join(beluga_main_folder, line[-2])
                            if line[-3]:
                                line[-3] = os.path.join(beluga_main_folder, line[-3])
                            readset_rna_out.append("\t".join(line))
                # Already analyzed sample RNA T
                to_transfer_rna_sample = rna_sample
                try:
                    analyzed_rna_sample = os.path.basename(glob.glob(os.path.join(beluga_main_raw_reads_folder, f"{patient}-*RT"))[0])
                except IndexError:
                    analyzed_rna_sample = None
                if analyzed_rna_sample:
                    for analysed_readset in glob.glob(os.path.join(beluga_main_raw_reads_folder, f"{patient}-*RT", "*_readset.tsv")):
                        with open(analysed_readset, 'rt') as readset_in:
                            reader = csv.reader(readset_in, delimiter="\t")
                            next(reader, None)
                            for line in reader:
                                # Needed to match topup to existing sample even if name is not the exact same
                                line[0] = analyzed_rna_sample
                                # Making fastq/bam paths absolute
                                if line[-1]:
                                    line[-1] = os.path.join(beluga_main_folder, "raw_reads", analyzed_rna_sample, os.path.basename(line[-1]))
                                if line[-2]:
                                    line[-2] = os.path.join(beluga_main_folder, "raw_reads", analyzed_rna_sample, os.path.basename(line[-2]))
                                if line[-3]:
                                    line[-3] = os.path.join(beluga_main_folder, "raw_reads", analyzed_rna_sample, os.path.basename(line[-3]))
                                readset_rna_out.append("\t".join(line))
                    rna_sample = analyzed_rna_sample
                if not args.dry_run:
                    # Move file
                    rna_files = os.listdir(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_rna_sample))
                    os.makedirs(os.path.join(beluga_main_raw_reads_folder, rna_sample), exist_ok=True)
                    for file_name in rna_files:
                        # print(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_rna_sample, file_name))
                        shutil.move(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_rna_sample, file_name), os.path.join(beluga_main_raw_reads_folder, rna_sample))
                    os.rmdir(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_rna_sample))
                    # Update database
                    Update_Samples_Table(
                        connection,
                        patient,
                        patient,
                        institution,
                        cohort,
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        rna_sample,
                        rna_sample
                        )
                # Writting outputs for each patient
                readset_rna_file = os.path.join(args.output_folder, "readset_files_RNA", f"{patient}_{date_formatted}_RNA_readset.tsv")
                os.makedirs(os.path.dirname(readset_rna_file), exist_ok=True)
                with open(readset_rna_file, "w", encoding="utf-8") as readset_file:
                    readset_file.write(f"{readset_header}\n")
                    for readset_line in readset_rna_out:
                        readset_file.write(f"{readset_line}\n")
                print(f"Generated {readset_rna_file}")
                readset_rna_out = []
                # else:
                #     print(" ".join([sample.sample for _, sample in sample.items()]))
            # if not args.dry_run:
            # with open(readset_rna_file, "w", encoding="utf-8") as readset_file:
            #     readset_file.write(f"{readset_header}\n")
            #     for readset_line in readset_rna_out:
            #         readset_file.write(f"{readset_line}\n")
            # for readset_line in readset_dna_out:
            #     print(readset_line)
            # print (f"Generated {readset_rna_file}")

        else:
            print("No RNA Samples to Move")

    # DNA
    elif sequencing_type == 'DNA':
        if patient_dict:
            readset_dna_out = []
            pair_out = []
            # if args.dry_run:
            #     print("The following sample(s) will be moved and a readset file + pair file will be created for a GenPipes analysis:")
            for patient, sample in patient_dict.items():
                sample_n = None
                sample_t = None
                try:
                    sample_n = sample['DN'].sample
                    patient = sample['DN'].patient
                    cohort = sample['DN'].cohort
                    institution = sample['DN'].institution
                except KeyError:
                    pass
                try:
                    sample_t = sample['DT'].sample
                    patient = sample['DT'].patient
                    cohort = sample['DT'].cohort
                    institution = sample['DT'].institution
                except KeyError:
                    pass
                # Pair
                if sample_n and sample_t:
                    transferred_readsets = glob.glob(os.path.join(beluga_transferred_raw_reads_folder, sample_n, "*_readset.tsv"))
                    if not transferred_readsets:
                        sys.exit(f"Sample {sample_n} doesn't have a readset file, there is an issue with the transfer. Please check {os.path.join(beluga_transferred_raw_reads_folder, sample_n)}\nExiting...")
                    for transferred_readset in transferred_readsets:
                        with open(transferred_readset, 'rt') as readset_in:
                            reader = csv.reader(readset_in, delimiter="\t")
                            next(reader, None)
                            for line in reader:
                                # Making fastq/bam paths absolute
                                if line[-1]:
                                    line[-1] = os.path.join(beluga_main_folder, line[-1])
                                if line[-2]:
                                    line[-2] = os.path.join(beluga_main_folder, line[-2])
                                if line[-3]:
                                    line[-3] = os.path.join(beluga_main_folder, line[-3])
                                readset_dna_out.append("\t".join(line))
                    transferred_readsets = glob.glob(os.path.join(beluga_transferred_raw_reads_folder, sample_t, "*_readset.tsv"))
                    if not transferred_readsets:
                        sys.exit(f"Sample {sample_t} doesn't have a readset file, there is an issue with the transfer. Please check {os.path.join(beluga_transferred_raw_reads_folder, sample_t)}\nExiting...")
                    for transferred_readset in transferred_readsets:
                        with open(transferred_readset, 'rt') as readset_in:
                            reader = csv.reader(readset_in, delimiter="\t")
                            next(reader, None)
                            for line in reader:
                                # Making fastq/bam paths absolute
                                if line[-1]:
                                    line[-1] = os.path.join(beluga_main_folder, line[-1])
                                if line[-2]:
                                    line[-2] = os.path.join(beluga_main_folder, line[-2])
                                if line[-3]:
                                    line[-3] = os.path.join(beluga_main_folder, line[-3])
                                readset_dna_out.append("\t".join(line))
                    # Already analyzed sample Normal
                    to_transfer_sample_n = sample_n
                    try:
                        analyzed_sample_n = os.path.basename(glob.glob(os.path.join(beluga_main_raw_reads_folder, f"{patient}-*DN"))[0])
                    except IndexError:
                        analyzed_sample_n = None
                    if analyzed_sample_n:
                        for analysed_readset in glob.glob(os.path.join(beluga_main_raw_reads_folder, f"{patient}-*DN", "*_readset.tsv")):
                            with open(analysed_readset, 'rt') as readset_in:
                                reader = csv.reader(readset_in, delimiter="\t")
                                next(reader, None)
                                for line in reader:
                                    # Needed to match topup to existing sample even if name is not the exact same
                                    line[0] = analyzed_sample_n
                                    # Making fastq/bam paths absolute
                                    if line[-1]:
                                        line[-1] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-1]))
                                    if line[-2]:
                                        line[-2] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-2]))
                                    if line[-3]:
                                        line[-3] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-3]))
                                    readset_dna_out.append("\t".join(line))
                        sample_n = analyzed_sample_n
                    # Already analyzed sample Tumour
                    to_transfer_sample_t = sample_t
                    try:
                        analyzed_sample_t = os.path.basename(glob.glob(os.path.join(beluga_main_raw_reads_folder, f"{patient}-*DT"))[0])
                    except IndexError:
                        analyzed_sample_t = None
                    if analyzed_sample_t:
                        for analysed_readset in glob.glob(os.path.join(beluga_main_raw_reads_folder, f"{patient}-*DT", "*_readset.tsv")):
                            with open(analysed_readset, 'rt') as readset_in:
                                reader = csv.reader(readset_in, delimiter="\t")
                                next(reader, None)
                                for line in reader:
                                    # Needed to match topup to existing sample even if name is not the exact same
                                    line[0] = analyzed_sample_t
                                    # Making fastq/bam paths absolute
                                    if line[-1]:
                                        line[-1] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-1]))
                                    if line[-2]:
                                        line[-2] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-2]))
                                    if line[-3]:
                                        line[-3] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-3]))
                                    readset_dna_out.append("\t".join(line))
                        sample_t = analyzed_sample_t
                    pair_out.append(f"{patient},{sample_n},{sample_t}")
                    if not args.dry_run:
                        # Move file
                        dna_n_files = os.listdir(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_sample_n))
                        os.makedirs(os.path.join(beluga_main_raw_reads_folder, sample_n), exist_ok=True)
                        for file_name in dna_n_files:
                            # print(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_sample_n, file_name))
                            shutil.move(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_sample_n, file_name), os.path.join(beluga_main_raw_reads_folder, sample_n))
                        os.rmdir(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_sample_n))
                        dna_t_files = os.listdir(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_sample_t))
                        os.makedirs(os.path.join(beluga_main_raw_reads_folder, sample_t), exist_ok=True)
                        for file_name in dna_t_files:
                            # print(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_sample_t, file_name))
                            shutil.move(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_sample_t, file_name), os.path.join(beluga_main_raw_reads_folder, sample_t))
                        os.rmdir(os.path.join(beluga_transferred_raw_reads_folder, to_transfer_sample_t))
                        # Update database
                        Update_Samples_Table(
                            connection,
                            patient,
                            patient,
                            institution,
                            cohort,
                            sample_n,
                            sample_n,
                            sample_t,
                            sample_t,
                            "NA",
                            "NA"
                            )
                    # else:
                    #     print(" ".join([sample.sample for _, sample in sample.items()]))
                # Topup not pair
                elif sample_n:
                    try:
                        analyzed_sample_n = os.path.basename(glob.glob(os.path.join(beluga_main_raw_reads_folder, f"{patient}-*DN"))[0])
                    except IndexError:
                        # Set to None as string to make the glob below not find without failing
                        analyzed_sample_n = "None"
                    if glob.glob(os.path.join(beluga_main_raw_reads_folder, analyzed_sample_n, "*.bam")):
                        # if not args.dry_run:
                        transferred_readsets = glob.glob(os.path.join(beluga_transferred_raw_reads_folder, sample_n, "*_readset.tsv"))
                        if not transferred_readsets:
                            sys.exit(f"Sample {sample_n} doesn't have a readset file, there is an issue with the transfer. Please check {os.path.join(beluga_transferred_raw_reads_folder, sample_n)}\nExiting...")
                        for transferred_readset in transferred_readsets:
                            with open(transferred_readset, 'rt') as readset_in:
                                reader = csv.reader(readset_in, delimiter="\t")
                                next(reader, None)
                                for line in reader:
                                    # Needed to match topup to existing sample even if name is not the exact same
                                    line[0] = analyzed_sample_n
                                    # Making fastq/bam paths absolute
                                    if line[-1]:
                                        line[-1] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-1]))
                                    if line[-2]:
                                        line[-2] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-2]))
                                    if line[-3]:
                                        line[-3] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-3]))
                                    readset_dna_out.append("\t".join(line))
                        for analysed_readset in glob.glob(os.path.join(beluga_main_raw_reads_folder, analyzed_sample_n, "*_readset.tsv")):
                            with open(analysed_readset, 'rt') as readset_in:
                                reader = csv.reader(readset_in, delimiter="\t")
                                next(reader, None)
                                for line in reader:
                                    # Needed to match topup to existing sample even if name is not the exact same
                                    line[0] = analyzed_sample_n
                                    # Making fastq/bam paths absolute
                                    if line[-1]:
                                        line[-1] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-1]))
                                    if line[-2]:
                                        line[-2] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-2]))
                                    if line[-3]:
                                        line[-3] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-3]))
                                    readset_dna_out.append("\t".join(line))
                        final_sample_n = analyzed_sample_n
                    else:
                        transferred_readsets = glob.glob(os.path.join(beluga_transferred_raw_reads_folder, sample_n, "*_readset.tsv"))
                        if not transferred_readsets:
                            sys.exit(f"Sample {sample_n} doesn't have a readset file, there is an issue with the transfer. Please check {os.path.join(beluga_transferred_raw_reads_folder, sample_n)}\nExiting...")
                        for transferred_readset in transferred_readsets:
                            with open(transferred_readset, 'rt') as readset_in:
                                reader = csv.reader(readset_in, delimiter="\t")
                                next(reader, None)
                                for line in reader:
                                    # Making fastq/bam paths absolute
                                    if line[-1]:
                                        line[-1] = os.path.join(beluga_main_folder, line[-1])
                                    if line[-2]:
                                        line[-2] = os.path.join(beluga_main_folder, line[-2])
                                    if line[-3]:
                                        line[-3] = os.path.join(beluga_main_folder, line[-3])
                                    readset_dna_out.append("\t".join(line))
                        final_sample_n = sample_n
                    # Already analyzed sample Tumor
                    try:
                        analyzed_sample_t = os.path.basename(glob.glob(os.path.join(beluga_main_raw_reads_folder, f"{patient}-*DT"))[0])
                        for analysed_readset in glob.glob(os.path.join(beluga_main_raw_reads_folder, analyzed_sample_t, "*_readset.tsv")):
                            with open(analysed_readset, 'rt') as readset_in:
                                reader = csv.reader(readset_in, delimiter="\t")
                                next(reader, None)
                                for line in reader:
                                    # Needed to match topup to existing sample even if name is not the exact same
                                    line[0] = analyzed_sample_t
                                    # Making fastq/bam paths absolute
                                    if line[-1]:
                                        line[-1] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-1]))
                                    if line[-2]:
                                        line[-2] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-2]))
                                    if line[-3]:
                                        line[-3] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-3]))
                                    readset_dna_out.append("\t".join(line))
                    except IndexError:
                        analyzed_sample_t = None
                    if analyzed_sample_t and final_sample_n:
                        pair_out.append(f"{patient},{analyzed_sample_n},{analyzed_sample_t}")
                        if not args.dry_run:
                            # Move file
                            dna_n_files = os.listdir(os.path.join(beluga_transferred_raw_reads_folder, sample_n))
                            os.makedirs(os.path.join(beluga_main_raw_reads_folder, analyzed_sample_n), exist_ok=True)
                            for file_name in dna_n_files:
                                # print(os.path.join(beluga_transferred_raw_reads_folder, sample_n, file_name))
                                shutil.move(os.path.join(beluga_transferred_raw_reads_folder, sample_n, file_name), os.path.join(beluga_main_raw_reads_folder, analyzed_sample_n))
                            os.rmdir(os.path.join(beluga_transferred_raw_reads_folder, sample_n))
                    # else:
                    #     print(" ".join([sample.sample for _, sample in sample.items()]))
                elif sample_t:
                    try:
                        analyzed_sample_t = os.path.basename(glob.glob(os.path.join(beluga_main_raw_reads_folder, f"{patient}-*DT"))[0])
                    except IndexError:
                        # Set to None as string to make the glob below not find without failing
                        analyzed_sample_t = "None"
                    if glob.glob(os.path.join(beluga_main_raw_reads_folder, analyzed_sample_t, "*.bam")):
                        # if not args.dry_run:
                        transferred_readsets = glob.glob(os.path.join(beluga_transferred_raw_reads_folder, sample_t, "*_readset.tsv"))
                        if not transferred_readsets:
                            sys.exit(f"Sample {sample_t} doesn't have a readset file, there is an issue with the transfer. Please check {os.path.join(beluga_transferred_raw_reads_folder, sample_t)}\nExiting...")
                        for transferred_readset in transferred_readsets:
                            with open(transferred_readset, 'rt') as readset_in:
                                reader = csv.reader(readset_in, delimiter="\t")
                                next(reader, None)
                                for line in reader:
                                    # Needed to match topup to existing sample even if name is not the exact same
                                    line[0] = analyzed_sample_t
                                    # Making fastq/bam paths absolute
                                    if line[-1]:
                                        line[-1] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-1]))
                                    if line[-2]:
                                        line[-2] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-2]))
                                    if line[-3]:
                                        line[-3] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-3]))
                                    readset_dna_out.append("\t".join(line))
                        for analysed_readset in glob.glob(os.path.join(beluga_main_raw_reads_folder, analyzed_sample_t, "*_readset.tsv")):
                            with open(analysed_readset, 'rt') as readset_in:
                                reader = csv.reader(readset_in, delimiter="\t")
                                next(reader, None)
                                for line in reader:
                                    # Needed to match topup to existing sample even if name is not the exact same
                                    line[0] = analyzed_sample_t
                                    # Making fastq/bam paths absolute
                                    if line[-1]:
                                        line[-1] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-1]))
                                    if line[-2]:
                                        line[-2] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-2]))
                                    if line[-3]:
                                        line[-3] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_t, os.path.basename(line[-3]))
                                    readset_dna_out.append("\t".join(line))
                        final_sample_t = analyzed_sample_t
                    else:
                        transferred_readsets = glob.glob(os.path.join(beluga_transferred_raw_reads_folder, sample_t, "*_readset.tsv"))
                        if not transferred_readsets:
                            sys.exit(f"Sample {sample_t} doesn't have a readset file, there is an issue with the transfer. Please check {os.path.join(beluga_transferred_raw_reads_folder, sample_t)}\nExiting...")
                        for transferred_readset in transferred_readsets:
                            with open(transferred_readset, 'rt') as readset_in:
                                reader = csv.reader(readset_in, delimiter="\t")
                                next(reader, None)
                                for line in reader:
                                    # Making fastq/bam paths absolute
                                    if line[-1]:
                                        line[-1] = os.path.join(beluga_main_folder, line[-1])
                                    if line[-2]:
                                        line[-2] = os.path.join(beluga_main_folder, line[-2])
                                    if line[-3]:
                                        line[-3] = os.path.join(beluga_main_folder, line[-3])
                                    readset_dna_out.append("\t".join(line))
                        final_sample_t = sample_t
                    # Already analyzed sample Normal
                    try:
                        analyzed_sample_n = os.path.basename(glob.glob(os.path.join(beluga_main_raw_reads_folder, f"{patient}-*DN"))[0])
                        for analysed_readset in glob.glob(os.path.join(beluga_main_raw_reads_folder, analyzed_sample_n, "*_readset.tsv")):
                            with open(analysed_readset, 'rt') as readset_in:
                                reader = csv.reader(readset_in, delimiter="\t")
                                next(reader, None)
                                for line in reader:
                                    # Needed to match topup to existing sample even if name is not the exact same
                                    line[0] = analyzed_sample_n
                                    # Making fastq/bam paths absolute
                                    if line[-1]:
                                        line[-1] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-1]))
                                    if line[-2]:
                                        line[-2] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-2]))
                                    if line[-3]:
                                        line[-3] = os.path.join(beluga_main_folder, "raw_reads", analyzed_sample_n, os.path.basename(line[-3]))
                                    readset_dna_out.append("\t".join(line))
                    except IndexError:
                        analyzed_sample_n = None
                    if analyzed_sample_n and final_sample_t:
                        pair_out.append(f"{patient},{analyzed_sample_n},{final_sample_t}")
                        if not args.dry_run:
                            # Move file
                            dna_t_files = os.listdir(os.path.join(beluga_transferred_raw_reads_folder, sample_t))
                            os.makedirs(os.path.join(beluga_main_raw_reads_folder, analyzed_sample_t), exist_ok=True)
                            for file_name in dna_t_files:
                                # print(os.path.join(beluga_transferred_raw_reads_folder, sample_t, file_name))
                                shutil.move(os.path.join(beluga_transferred_raw_reads_folder, sample_t, file_name), os.path.join(beluga_main_raw_reads_folder, analyzed_sample_t))
                            os.rmdir(os.path.join(beluga_transferred_raw_reads_folder, sample_t))
                # Writting outputs for each patient
                if readset_dna_out and pair_out:
                    readset_dna_file = os.path.join(args.output_folder, "readset_pair_files_DNA", f"{patient}_{date_formatted}_TP_readset.tsv")
                    pair_dna_file = os.path.join(args.output_folder, "readset_pair_files_DNA", f"{patient}_{date_formatted}_TP_pairs.csv")
                    os.makedirs(os.path.dirname(readset_dna_file), exist_ok=True)
                    os.makedirs(os.path.dirname(pair_dna_file), exist_ok=True)
                    with open(readset_dna_file, "w", encoding="utf-8") as readset_file, open(pair_dna_file, "w", encoding="utf-8") as pair_file:
                        readset_file.write(f"{readset_header}\n")
                        for readset_line in readset_dna_out:
                            readset_file.write(f"{readset_line}\n")
                        for pair_line in pair_out:
                            pair_file.write(f"{pair_line}\n")
                    print(f"Generated {readset_dna_file} and {pair_dna_file}")
                    readset_dna_out = []
                    pair_out = []

            # if not args.dry_run:
            # with open(readset_dna_file, "w", encoding="utf-8") as readset_file, open(pair_dna_file, "w", encoding="utf-8") as pair_file:
            #     readset_file.write(f"{readset_header}\n")
            #     for readset_line in readset_dna_out:
            #         readset_file.write(f"{readset_line}\n")
            #     for pair_line in pair_out:
            #         pair_file.write(f"{pair_line}\n")
            # for readset_line in readset_dna_out:
            #     print(readset_line)
            # for pair_line in pair_out:
            #     print(pair_line)
            # print(f"Generated {readset_dna_file} and {pair_dna_file}")

        else:
            print("No DNA pairs to Move")


    connection.commit()
    connection.close()

class Sample():
    """docstring for Sample"""
    def __init__(self, sample):
        super(Sample, self).__init__()
        self.sample = sample
        result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
        self.patient = result.group(1)
        self.cohort = result.group(2)
        self.institution = result.group(3)
        self.type = result.group(4) + result.group(5)

def intersection(lst1, lst2):
    """Intersection of 2 lists"""
    # Use of hybrid method
    temp = set(lst2)
    return [value for value in lst1 if value in temp]


if __name__ == '__main__':
    main()
