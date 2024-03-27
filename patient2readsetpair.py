#!/usr/bin/env python3

import argparse
import glob
import csv
import os
import logging

logging.basicConfig(format='%(levelname)s: %(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(prog='tmp.py', description="Creates readset and pair files for DNA")
    parser.add_argument('--input', required=True, help="Input file with list of patients")
    parser.add_argument('--output', required=False, help="Output prefix (Default: <input_filename>)")
    args = parser.parse_args()

    if not args.output:
        output = f"{os.path.basename(args.input).split('.')[0]}"
    else:
        output = args.output

    readsetpair(args.input, output)

def readsetpair(input, output):
    # beluga_moh_folder = "/Users/pstretenowich/Mount_points/beluga/MOH_PROCESSING"
    beluga_moh_folder = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING"
    beluga_main_folder = os.path.join(beluga_moh_folder, "MAIN")
    beluga_rr_folder = os.path.join(beluga_main_folder, "raw_reads")

    readset_header = 'Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM'

    pair_out = []
    readset_dna_out = []

    with open(input, 'r') as patients:
        for patient in patients:
            patient = patient.strip()
            sample_t = glob.glob(os.path.join(beluga_rr_folder, f"{patient}-*DT"))[0].split("/")[-1]
            sample_n = glob.glob(os.path.join(beluga_rr_folder, f"{patient}-*DN"))[0].split("/")[-1]
            for analysed_readset in glob.glob(os.path.join(beluga_rr_folder, sample_t, "*_readset.tsv")):
                with open(analysed_readset, 'rt') as readset_in:
                    reader = csv.reader(readset_in, delimiter="\t")
                    next(reader, None)
                    for line in reader:
                        # Needed to match topup to existing sample even if name is not the exact same
                        line[0] = sample_t
                        # Making fastq/bam paths absolute
                        if line[-1]:
                            line[-1] = os.path.join(beluga_main_folder, line[-1])
                        if line[-2]:
                            line[-2] = os.path.join(beluga_main_folder, line[-2])
                        if line[-3]:
                            line[-3] = os.path.join(beluga_main_folder, line[-3])
                        readset_dna_out.append("\t".join(line))
            for analysed_readset in glob.glob(os.path.join(beluga_rr_folder, sample_n, "*_readset.tsv")):
                with open(analysed_readset, 'rt') as readset_in:
                    reader = csv.reader(readset_in, delimiter="\t")
                    next(reader, None)
                    for line in reader:
                        # Needed to match topup to existing sample even if name is not the exact same
                        line[0] = sample_n
                        # Making fastq/bam paths absolute
                        if line[-1]:
                            line[-1] = os.path.join(beluga_main_folder, line[-1])
                        if line[-2]:
                            line[-2] = os.path.join(beluga_main_folder, line[-2])
                        if line[-3]:
                            line[-3] = os.path.join(beluga_main_folder, line[-3])
                        readset_dna_out.append("\t".join(line))
            pair_out.append(f"{patient},{sample_n},{sample_t}")
            # glob.glob(os.path.join(beluga_transferred_raw_reads_folder, sample_n, "*_readset.tsv"))

    readset_dna_file = f"{output}_TP_readset.tsv"
    pair_dna_file = f"{output}_TP_pairs.csv"
    with open(readset_dna_file, "w", encoding="utf-8") as readset_file, open(pair_dna_file, "w", encoding="utf-8") as pair_file:
        readset_file.write(f"{readset_header}\n")
        for readset_line in readset_dna_out:
            readset_file.write(f"{readset_line}\n")
        for pair_line in pair_out:
            pair_file.write(f"{pair_line}\n")



if __name__ == '__main__':
    main()
