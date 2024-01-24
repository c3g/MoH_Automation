#!/usr/bin/env python3

import glob
import os
import re
import csv
import progressbar

WIDGETS = [' [', progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') - ', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']

class Sample():
    """docstring for Sample"""
    def __init__(self, sample):
        super(Sample, self).__init__()
        self.sample = sample
        result = re.search(r"^((MoHQ-(JG|HM|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
        self.patient = result.group(1)
        self.cohort = result.group(2)
        self.institution = result.group(3)
        self.type = result.group(4) + result.group(5)

def main():
    main_beluga_folder = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN"
    main_raw_reads_folder = os.path.join(main_beluga_folder, "raw_reads")
    readset_header = 'Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM'

    # Populate the lists
    patient_dict = {}
    dna_patient_pair_dict = {}
    all_samples_list = os.listdir(main_raw_reads_folder)
    print("Finding all samples...")
    with progressbar.ProgressBar(max_value=len(all_samples_list), widgets=WIDGETS) as progress:
        for index, sample in enumerate(all_samples_list, 1):
        # for sample in all_samples_list:
            if sample.startswith("MoHQ"):
                sample_object = Sample(sample)
                if sample.endswith('RT'):
                    try:
                        patient_dict[sample_object.patient][sample_object.type] = sample_object
                        # dna_patient_pair_dict[sample_object.patient] = patient_dict[sample_object.patient]
                    except KeyError:
                        patient_dict[sample_object.patient] = {sample_object.type: sample_object}
                elif sample.endswith('DT') or sample.endswith('DN'):
                    try:
                        patient_dict[sample_object.patient][sample_object.type] = sample_object
                        dna_patient_pair_dict[sample_object.patient] = patient_dict[sample_object.patient]
                    except KeyError:
                        patient_dict[sample_object.patient] = {sample_object.type: sample_object}
            progress.update(index)

    readset_rna_out = []
    readset_rna_file = os.path.join(main_beluga_folder, "all_RNA_readset.tsv")
    readset_dna_out = []
    pair_out = []
    readset_dna_file = os.path.join(main_beluga_folder, "all_TP_readset.tsv")
    pair_dna_file = os.path.join(main_beluga_folder, "all_TP_pairs.csv")
    print("Gathering readsets/pair for all samples per sequencing type (TP/RNA)...")
    index = 1
    with progressbar.ProgressBar(max_value=len(patient_dict), widgets=WIDGETS) as progress:
        for patient, sample in patient_dict.items():
            sample_n = None
            sample_t = None
            sample_rna = None
            try:
                sample_n = sample['DN'].sample
                patient = sample['DN'].patient
            except KeyError:
                pass
            try:
                sample_t = sample['DT'].sample
                patient = sample['DT'].patient
            except KeyError:
                pass
            try:
                sample_rna = sample['RT'].sample
                patient = sample['RT'].patient
            except KeyError:
                pass

            if sample_rna:
                sample_readsets = glob.glob(os.path.join(main_raw_reads_folder, sample_rna, "*_readset.tsv"))
                for sample_readset in sample_readsets:
                    with open(sample_readset, 'rt') as readset_in:
                        reader = csv.reader(readset_in, delimiter="\t")
                        next(reader, None)
                        for line in reader:
                            readset_rna_out.append("\t".join(line))

            if sample_n and sample_t:
                sample_readsets = glob.glob(os.path.join(main_raw_reads_folder, sample_n, "*_readset.tsv"))
                for sample_readset in sample_readsets:
                    with open(sample_readset, 'rt') as readset_in:
                        reader = csv.reader(readset_in, delimiter="\t")
                        next(reader, None)
                        for line in reader:
                            readset_dna_out.append("\t".join(line))
                sample_readsets = glob.glob(os.path.join(main_raw_reads_folder, sample_t, "*_readset.tsv"))
                for sample_readset in sample_readsets:
                    with open(sample_readset, 'rt') as readset_in:
                        reader = csv.reader(readset_in, delimiter="\t")
                        next(reader, None)
                        for line in reader:
                            readset_dna_out.append("\t".join(line))
                pair_out.append(f"{patient},{sample_n},{sample_t}")

            progress.update(index)
            index += 1

    print("Writing readsets/pair for all samples per sequencing type (TP/RNA)...")
    # For RNA
    with open(readset_rna_file, "w", encoding="utf-8") as readset_file:
        readset_file.write(f"{readset_header}\n")
        for readset_line in readset_rna_out:
            readset_file.write(f"{readset_line}\n")

    # For DNA
    with open(readset_dna_file, "w", encoding="utf-8") as readset_file, open(pair_dna_file, "w", encoding="utf-8") as pair_file:
        readset_file.write(f"{readset_header}\n")
        for readset_line in readset_dna_out:
            readset_file.write(f"{readset_line}\n")
        for pair_line in pair_out:
            pair_file.write(f"{pair_line}\n")

    print("...Done!")

    # dna_patient_pair_dict = {}
    # # rna_patient_dict = {}
    # dna_patient_dict = {}
    # for sample in glob.glob(os.path.join(main_raw_reads_folder, "*")):
    #     sample = sample.split("/")[-1]
    #     if sample != "bad":
    #         sample_object = Sample(sample)
    #         sample_readsets = glob.glob(os.path.join(main_raw_reads_folder, sample, "*_readset.tsv"))
    #         for sample_readset in sample_readsets:
    #             with open(sample_readset, 'rt') as readset_in:
    #                 reader = csv.reader(readset_in, delimiter="\t")
    #                 next(reader, None)
    #                 for line in reader:
    #                     if sample.endswith('RT'):
    #                         rna_line.append("\t".join(line))
    #                     elif sample.endswith('DT') or sample.endswith('DN'):
    #                         dna_line.append("\t".join(line))
    #                         try:
    #                             dna_patient_dict[sample_object.patient][sample_object.type] = readset_line
    #                             dna_patient_pair_dict[sample_object.patient] = dna_patient_dict[sample_object.patient]
    #                         except KeyError:
    #                             dna_patient_dict[sample_object.patient] = {sample_object.type: readset_line}

            # print(sample)
            # sample_object = Sample(sample)
            # sample_readset = os.path.join(main_raw_reads_folder, sample, f"{sample}_readset.tsv")
            # with open(sample_readset, "r", encoding="utf-8") as file:
            #     readset_line = file.readlines()[1:]
            #     readset_line = "".join(readset_line).strip()
            #     if sample.endswith('RT'):
            #         # print(sample)
            #         rna_line.append(readset_line)
            #     elif sample.endswith('DT') or sample.endswith('DN'):
            #         dna_line.append(readset_line.split("\t")[0])
            #         try:
            #             dna_patient_dict[sample_object.patient][sample_object.type] = readset_line
            #             dna_patient_pair_dict[sample_object.patient] = dna_patient_dict[sample_object.patient]
            #             # print(sample)
            #             # dna_line.append(readset_line)
            #         except KeyError:
            #             dna_patient_dict[sample_object.patient] = {sample_object.type: readset_line}
            #             # print(sample)
            #             # dna_line.append(readset_line)



    # for patient, sample in dna_patient_pair_dict.items():
    #     for dna_sample in dna_line:
    #         if sample['DN'].split("\t")[0] == dna_sample:
    #             dna_line.remove(dna_sample)
    #         if sample['DT'].split("\t")[0] == dna_sample:
    #             dna_line.remove(dna_sample)
    # print(dna_line)


    # with open(os.path.join(main_beluga_folder, "all_RNA_readset.tsv"), "w", encoding="utf-8") as file:
    #     file.write(f"{readset_header}\n")
    #     for line in rna_line:
    #         file.write(f"{line}\n")

    # with open(os.path.join(main_beluga_folder, "all_TP_readset.tsv"), "w", encoding="utf-8") as readset_file, open(os.path.join(main_beluga_folder, "all_TP_pairs.csv"), "w", encoding="utf-8") as pair_file:
    #     readset_file.write(f"{readset_header}\n")
    #     # for line in dna_line:
    #     #     readset_file.write(f"{line}\n")
    #     for patient, sample in dna_patient_pair_dict.items():
    #         sample_n = sample['DN'].split("\t")[0]
    #         sample_t = sample['DT'].split("\t")[0]
    #         # sample_n = sample['DN'].sample
    #         # sample_t = sample['DT'].sample
    #         readset_file.write(f"{sample['DN']}\n")
    #         readset_file.write(f"{sample['DT']}\n")
    #         pair_file.write(f"{patient},{sample_n},{sample_t}\n")


if __name__ == '__main__':
    main()
