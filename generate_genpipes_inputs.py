#!/usr/bin/env python3

import glob
import sys
import os
import re
import datetime
import shutil
import argparse
from  moh_resources import create_connection, update_patient_table, extract_readset_details

def main():
    parser = argparse.ArgumentParser(prog='generate_genpipes_inputs.py', description="Generate readset file for GenPipes and move raw_reads in right folder for MOH project.")
    parser.add_argument('--type', required=True, help="Type of analysis: either RNA or DNA", choices=['DNA', 'RNA'])
    parser.add_argument('--dry_run', required=False, help="Will not move raw_data reads but prints GenPipes readset", default=False, action='store_true')
    parser.add_argument('--force', required=False, help="If used the duplicates samples will be overwritten", default=False, action='store_true')
    args = parser.parse_args()

    sequencing_type = args.type

    beluga_moh_folder = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING"
    beluga_main_folder = os.path.join(beluga_moh_folder, "MAIN")
    # beluga_db = os.path.join(beluga_moh_folder, "DATABASE", MOH_analysis.db)
    # TEST
    beluga_db = "/Users/pstretenowich/Documents/local/projects/MOH/MOH_analysis.db"
    beluga_main_raw_reads_folder = os.path.join(beluga_main_folder, "raw_reads")
    beluga_transferred_raw_reads_folder = os.path.join(beluga_moh_folder, "raw_reads")

    date = datetime.datetime.today()
    date_formatted = date.strftime("%Y-%m-%dT%H.%M.%S")
    #Connect to the db
    connection = create_connection(beluga_db)

    readset_header = 'Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM'

    # TODO: test writing to files, with /n and without

    # Populate the lists
    rna_samples_list = []
    patient_dict = {}
    dna_patient_pair_dict = {}
    transferred_samples_list = os.listdir(beluga_transferred_raw_reads_folder)
    duplicates = []
    for sample in transferred_samples_list:
        sample_object = Sample(sample)
        if sequencing_type == 'RNA' and sample.endswith('RT'):
            if os.path.exists(os.path.join(beluga_main_raw_reads_folder, sample)):
                duplicates.append(sample)
            # rna_samples_list.append(sample)
            try:
                patient_dict[sample_object.patient][sample_object.type] = sample_object
                # dna_patient_pair_dict[sample_object.patient] = patient_dict[sample_object.patient]
            except KeyError:
                patient_dict[sample_object.patient] = {sample_object.type: sample_object}
        elif sequencing_type == 'DNA' and (sample.endswith('DT') or sample.endswith('DN')):
            if os.path.exists(os.path.join(beluga_main_raw_reads_folder, sample)):
                duplicates.append(sample)
            # sample_object = Sample(sample)
            try:
                patient_dict[sample_object.patient][sample_object.type] = sample_object
                dna_patient_pair_dict[sample_object.patient] = patient_dict[sample_object.patient]
            except KeyError:
                patient_dict[sample_object.patient] = {sample_object.type: sample_object}

    if duplicates and not args.force:
        duplicates_file = os.path.join(beluga_main_folder, f"{date_formatted}_duplicates.txt")
        with open(duplicates_file, "w+", encoding="utf-8") as filename:
            for line in duplicates:
                filename.write(f"{line}\n")
        sys.exit(f"Samples are already present in in MAIN/raw_reads folder. Either use --force option to overwritte them or manually delete them. Cf. file {duplicates_file} for a liust of duplicates.\nExiting...")

    #RNA
    if sequencing_type == 'RNA':
        if patient_dict:
            # Check to see if any files are already present in the final directory
            # duplicate_samples(beluga_main_folder, beluga_main_raw_reads_folder, rna_samples_list, date_formatted)
            # readset generation
            readset_rna_file = os.path.join(beluga_main_folder, date_formatted + "_RNA_readset.tsv")
            # readset = []
            # for sample in rna_samples_list:
                # mylines = []
                # with open(os.path.join(beluga_transferred_raw_reads_folder, sample, sample + '_readset.tsv'), 'rt', encoding="utf-8") as myfile:
                #     for myline in myfile:
                #         mylines.append(myline)
                # readset = readset + mylines[1:]
            if args.dry_run:
                print("The following samples will be moved and a readset file will be create for a GenPipes analysis:")
            else:
                with open(readset_rna_file, "w", encoding="utf-8") as readset_file:
                    readset_file.write(readset_header)
            with open(readset_rna_file, "a+", encoding="utf-8") as readset_file:
                for patient, sample in patient_dict.items():
                    if not args.dry_run:
                        rna_sample = sample.sample
                        cohort = sample.cohort
                        institution = sample.institution
                        # Writting readset
                        readset_line = extract_readset_details(connection, rna_sample)
                        readset_line = [str(item or '') for item in readset_line]
                        readset_file.write("\t".join(readset_line))
                        # Updating db
                        update_patient_table(
                            connection,
                            patient,
                            patient,
                            institution,
                            cohort,
                            None,
                            None,
                            None,
                            None,
                            rna_sample,
                            rna_sample,
                            )
                        # move file
                        shutil.move(os.path.join(beluga_transferred_raw_reads_folder, sample), os.path.join(beluga_main_raw_reads_folder, sample))
                    else:
                        print(sample)
            if not args.dry_run:
                print (f"Generated {readset_rna_file}")

            #move files
            # for sample in rna_samples_list:
            #     shutil.move(os.path.join(beluga_transferred_raw_reads_folder, sample), os.path.join(beluga_main_raw_reads_folder, sample))

            #Add To database
            # for sample in rna_samples_list:
            #     result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
            #     patient = result.group(1)
            #     cohort = result.group(2)
            #     institution = result.group(3)
            #     rna_sample = sample
            #     update_patient_table(
            #         connection,
            #         patient,
            #         patient,
            #         institution,
            #         cohort,
            #         None,
            #         None,
            #         None,
            #         None,
            #         rna_sample,
            #         rna_sample,
            #         )

        else:
            print("No RNA Samples to Move")

    #DNA
    elif sequencing_type == 'DNA':
        if dna_patient_pair_dict:
            readset_dna_file = os.path.join(beluga_main_folder, date_formatted + "_TP_readset.tsv")
            pair_dna_file = os.path.join(beluga_main_folder, date_formatted + "_TP_pairs.csv")
            if args.dry_run:
                print("The following samples will be moved and a readset/pair file will be create for a GenPipes analysis:")
            else:
                with open(readset_dna_file, "w", encoding="utf-8") as readset_file:
                    readset_file.write(readset_header)
            with open(readset_dna_file, "a+", encoding="utf-8") as readset_file, open(pair_dna_file, "a+", encoding="utf-8") as pair_file:
                for patient, sample in dna_patient_pair_dict.items():
                    if not args.dry_run:
                        sample_n = sample['DN'].sample
                        sample_t = sample['DT'].sample
                        cohort = sample['DT'].cohort
                        institution = sample['DT'].institution
                        # Writting readset
                        readset_line_n = extract_readset_details(connection, sample_n)
                        readset_line_n = [str(item or '') for item in readset_line_n]
                        readset_file.write("\t".join(readset_line_n))
                        readset_line_t = extract_readset_details(connection, sample_t)
                        readset_line_t = [str(item or '') for item in readset_line_t]
                        readset_file.write("\t".join(readset_line_t))
                        # Writting pair
                        pair_file.write(f"{patient},{sample_n},{sample_t}")
                        # Updating db
                        update_patient_table(
                            connection,
                            patient,
                            patient,
                            institution,
                            cohort,
                            sample_n,
                            sample_n,
                            sample_t,
                            sample_t,
                            None,
                            None,
                            )
                        # move file
                        shutil.move(os.path.join(beluga_transferred_raw_reads_folder, sample_n), os.path.join(beluga_main_raw_reads_folder, sample_n))
                        shutil.move(os.path.join(beluga_transferred_raw_reads_folder, sample_t), os.path.join(beluga_main_raw_reads_folder, sample_t))
                    else:
                        print(sample)
            if not args.dry_run:
                print (f"Generated {readset_dna_file} and {pair_dna_file}")

            # # Check to see if any files are already present in the final directory
            # # duplicate_samples(beluga_main_folder, beluga_main_raw_reads_folder, dna_samples_list, date_formatted)
            # readset = []
            # readset.append('Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM\n')
            # pairs = []
            # for _, dna_pair in dna_patient_pair_dict.items():
            #     sample_n = dna_pair['DN'].sample
            #     sample_t = dna_pair['DT'].sample
            #     #readset
            #     readset_lines = []
            #     with open(os.path.join(beluga_transferred_raw_reads_folder, sample_n, sample_n + '_readset.tsv'), 'rt', encoding="utf-8") as file:
            #         for line in file:
            #             readset_lines.append(line)
            #     readset = readset + readset_lines[1:]
            #     readset_lines = []
            #     with open(os.path.join(beluga_transferred_raw_reads_folder, sample_t, sample_t + '_readset.tsv'), 'rt', encoding="utf-8") as file:
            #         for line in file:
            #             readset_lines.append(line)
            #     readset = readset + readset_lines[1:]
            #    #Add to db
            #     result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample_n)
            #     patient = result.group(1)
            #     cohort = result.group(2)
            #     institution = result.group(3)
            #     # path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*'
            #     # for filename in glob.glob(path):
            #     #     with open(filename, 'r', encoding="utf-8") as file:
            #     #         for line in file:
            #     #             if patient in line:
            #     #                 fields = line.split(",")
            #     #                 run = fields[0].split("/")[7]
            #     update_patient_table(
            #         connection,
            #         patient,
            #         patient,
            #         institution,
            #         cohort,
            #         sample_n,
            #         sample_n,
            #         sample_t,
            #         sample_t,
            #         None,
            #         None,
            #         )

            #     #pairs
            #     pairs.append(f"{patient},{sample_n},{sample_t}")

            #     #move file
            #     shutil.move(os.path.join(beluga_transferred_raw_reads_folder, sample_n), os.path.join(beluga_main_raw_reads_folder, sample_n))
            #     shutil.move(os.path.join(beluga_transferred_raw_reads_folder, sample_t), os.path.join(beluga_main_raw_reads_folder, sample_t))

            # #write the files
            # with open(os.path.join(beluga_main_folder, date_formatted + "_TP_readset.tsv"), "w+", encoding="utf-8") as filename:
            #     for line in readset:
            #         filename.write(f"{line}")
            #     print (f"Generated {os.path.join(beluga_main_folder, date_formatted)}_TP_readset.tsv")
            # with open(os.path.join(beluga_main_folder, date_formatted + "_TP_pairs.csv"), "w+", encoding="utf-8") as filename:
            #     for line in pairs:
            #         filename.write(f"{line}\n")
            #     print (f"Generated {os.path.join(beluga_main_folder, date_formatted)}_TP_pairs.csv")

        else:
            print("No DNA pairs to Move")


    connection.commit()
    connection.close()

def duplicate_samples(beluga_main_folder, beluga_main_raw_reads_folder, samples_list, date_formatted):
    duplicates = []
    for sample in samples_list:
        if os.path.exists(os.path.join(beluga_main_raw_reads_folder, sample)):
            duplicates.append(sample)
    duplicates_file = os.path.join(beluga_main_folder, f"{date_formatted}_duplicates.txt")
    if duplicates:
        with open(duplicates_file, "w+", encoding="utf-8") as filename:
            for line in duplicates:
                filename.write(f"{line}\n")
        sys.exit(f"Samples are already present in in MAIN/raw_reads folder. Either use --force option to overwritte them or manually delete them. Cf. file {duplicates_file} for a liust of duplicates.\nExiting...")

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


if __name__ == '__main__':
    main()
