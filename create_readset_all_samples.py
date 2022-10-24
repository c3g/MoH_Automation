#!/usr/bin/env python3

import glob
import os
import re


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

def main():
    main_beluga_folder = "/Users/pstretenowich/Mount_points/beluga/C3G/projects/MOH_PROCESSING/MAIN"
    main_raw_reads_folder = os.path.join(main_beluga_folder, "raw_reads")
    readset_header = 'Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM'
    rna_line = []
    dna_line = []

    dna_patient_pair_dict = {}
    rna_patient_dict = {}
    dna_patient_dict = {}
    for sample in glob.glob(os.path.join(main_raw_reads_folder, "*")):
        sample = sample.split("/")[-1]
        if sample != "bad":
            # print(sample)
            sample_object = Sample(sample)
            sample_readset = os.path.join(main_raw_reads_folder, sample, f"{sample}_readset.tsv")
            with open(sample_readset, "r", encoding="utf-8") as file:
                readset_line = file.readlines()[1:]
                readset_line = "".join(readset_line).strip()
                if sample.endswith('RT'):
                    # print(sample)
                    rna_line.append(readset_line)
                elif sample.endswith('DT') or sample.endswith('DN'):
                    dna_line.append(readset_line.split("\t")[0])
                    try:
                        dna_patient_dict[sample_object.patient][sample_object.type] = readset_line
                        dna_patient_pair_dict[sample_object.patient] = dna_patient_dict[sample_object.patient]
                        # print(sample)
                        # dna_line.append(readset_line)
                    except KeyError:
                        dna_patient_dict[sample_object.patient] = {sample_object.type: readset_line}
                        # print(sample)
                        # dna_line.append(readset_line)



    # for patient, sample in dna_patient_pair_dict.items():
    #     for dna_sample in dna_line:
    #         if sample['DN'].split("\t")[0] == dna_sample:
    #             dna_line.remove(dna_sample)
    #         if sample['DT'].split("\t")[0] == dna_sample:
    #             dna_line.remove(dna_sample)
    # print(dna_line)


    with open(os.path.join(main_beluga_folder, "all_RNA_readset.tsv"), "w", encoding="utf-8") as file:
        file.write(f"{readset_header}\n")
        for line in rna_line:
            file.write(f"{line}\n")

    with open(os.path.join(main_beluga_folder, "all_TP_readset.tsv"), "w", encoding="utf-8") as readset_file, open(os.path.join(main_beluga_folder, "all_TP_pairs.csv"), "w", encoding="utf-8") as pair_file:
        readset_file.write(f"{readset_header}\n")
        # for line in dna_line:
        #     readset_file.write(f"{line}\n")
        for patient, sample in dna_patient_pair_dict.items():
            sample_n = sample['DN'].split("\t")[0]
            sample_t = sample['DT'].split("\t")[0]
            # sample_n = sample['DN'].sample
            # sample_t = sample['DT'].sample
            readset_file.write(f"{sample['DN']}\n")
            readset_file.write(f"{sample['DT']}\n")
            pair_file.write(f"{patient},{sample_n},{sample_t}\n")


if __name__ == '__main__':
    main()

