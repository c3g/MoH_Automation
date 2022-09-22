#!/usr/bin/env python3

import argparse
import subprocess
import datetime
import glob
import os
import re
import progressbar
from moh_resources import create_connection, update_readset_table, PasswordPromptAction, update_key_metric_table_run_metrics, raw_mean_coverage_check, raw_reads_count_check, raw_duplication_rate_check, extract_value, update_sample_table, update_fileloc_details_run_processing

WIDGETS = [' [', progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') - ', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']

def main():
    parser = argparse.ArgumentParser(prog='transfer_raw_data.py', description="Transfers bams/fastqs from Abacus to Beluga for MOH project.")
    parser.add_argument('--run_csv', required=True, help="path/to the csv file of the run you want to transfer")
    parser.add_argument('--user', required=True, help="The username to be used for the transfer")
    parser.add_argument('--password', dest='password', action=PasswordPromptAction, type=str, required=False, help="Password for connecting to Beluga via ssh")
    # parser.add_argument('--force', required=False, help="If used the transfer will be submitted even if the run has already been transferred")
    parser.add_argument('--dry_run', required=False, help="Will not submit the globus transfer but prints files to be transferred", default=False, action='store_true')
    args = parser.parse_args()

    # Beluga locations
    beluga_moh_folder = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING"
    beluga_raw_reads_dir = os.path.join(beluga_moh_folder, "raw_reads")
    # beluga_db = os.path.join(beluga_moh_folder, "DATABASE", MOH_analysis.db)
    # TEST
    beluga_db = "/scratch/stretenp/moh_test/MOH_analysis.db"

    # Abacus locations
    abacus_moh_folder = "/lb/project/mugqic/projects/MOH"
    # abacus_db = os.path.join(abacus_moh_folder, "DATABASE", "MOH_analysis.db")
    # TEST
    abacus_db = "/Users/pstretenowich/Documents/local/projects/MOH/MOH_analysis.db"
    # abacus_temp = "/lb/project/mugqic/projects/MOH/TEMP"
    # TEST
    abacus_temp = "/Users/pstretenowich/Documents/local/projects/MOH/TEMP"
    abacus_temp = os.path.join(abacus_moh_folder, "TEMP")
    timestamp = datetime.datetime.now()
    timestamp_formatted = timestamp.strftime("%Y-%m-%dT%H.%M.%S")

    # TODO: add locking database https://stackoverflow.com/questions/14272715/how-can-i-lock-a-sqlite-database on Beluga while being modified on Abacus
    # Rsyncing Beluga database:
    # subprocess.run(f"""rsync -u {args.user}@beluga.calculcanada.ca:{beluga_db} {abacus_db}""", check=True, shell=True)

    connection = create_connection(abacus_db)

    #Beluga Endpoint
    beluga_endpoint = '278b9bfe-24da-11e9-9fa2-0a06afd4a22e'
    #Narval Endpoint #not used
    narval_endpoint = 'a1713da6-098f-40e6-b3aa-034efe8b6e5b'
    #abacus Endpoint
    abacus_endpoint = '6c66d53d-a79d-11e8-96fa-0a6d4e044368'

    # Checking number of lines in run csv
    completed_process = subprocess.run(["wc", "-l", args.run_csv], capture_output=True, check=True, encoding="utf-8")
    # print(int(completed_process.stdout.strip().split(" ")[0])-1)

    # TODO: Add top-up capacity by reading sample name and matching patient name

    run_csv_nb_lines = int(completed_process.stdout.strip().split(" ")[0]) - 1
    empty_transfer_list_file = False
    print("Parsing run files...")
    with open(args.run_csv, "r", encoding="utf-8") as filename, open(args.run_csv[:-3] + "align_bwa_mem.csv", "r", encoding="utf-8") as align_bwa_mem:
        with progressbar.ProgressBar(max_value=run_csv_nb_lines, widgets=WIDGETS) as progress:
            # skipping header line
            for _ in range(1):
                next(filename)
            # Gathering readset info
            for index, line in enumerate(filename, 1):
                line_splitted = line.split(",")
                # Grabbing all info for database
                run_name = line_splitted[0]
                sample = line_splitted[6]
                run_type = line_splitted[3]
                run = line_splitted[1]
                lane = line_splitted[2]
                adapter1 = line_splitted[27]
                adapter2 = line_splitted[28].strip()
                quality_offset = "33"
                bed = None
                library_id = line_splitted[9]
                library_type = f"{run}_{lane}"
                readset = f"{sample}.{library_type}"
                for align_bwa_mem_line in align_bwa_mem:
                    align_bwa_mem_line_splitted = align_bwa_mem_line.split(",")
                    if align_bwa_mem_line_splitted[6] == sample:
                        raw_reads_count = align_bwa_mem_line_splitted[12]
                        raw_duplication_rate = align_bwa_mem_line_splitted[15]
                        raw_median_insert_size = align_bwa_mem_line_splitted[37]
                        raw_mean_insert_size = align_bwa_mem_line_splitted[38]
                        raw_mean_coverage = align_bwa_mem_line_splitted[41]

                tester = re.match(r"^MoHQ-(JG|CM|GC|MU|MR|XX)-\w+-\w+-\w+-\w+(D|R)(T|N)", sample)
                transfer_list_file = os.path.join(abacus_temp, f"{timestamp_formatted}_{run_name}_transfer_list.txt")
                # Emptying transfer_list_file
                if not empty_transfer_list_file and not args.dry_run:
                    with open(transfer_list_file, 'w', encoding="utf-8"):
                        empty_transfer_list_file = True
                if not tester:
                    print(f"Warning: {sample} is in inproper format and has to be renamed!")
                elif sample.endswith('DN') or sample.endswith('DT'):
                    abacus_bam = os.path.join(os.path.dirname(args.run_csv), f"Aligned.{lane}", "alignment", sample, f"run{library_type}", f"{sample}_{library_id}.sorted.bam")
                    readset_bam = os.path.join("raw_reads", sample, f"{sample}.bam")
                    fastq1 = None
                    fastq2 = None
                    if not os.path.isfile(abacus_bam):
                        print(f"Warning: BAM file {abacus_bam} doesn't exist for sample {sample}.")
                    else:
                        # if os.path.getsize(transfer_list_file) != 0 and not args.force:
                        #     raise FileExistsError(f"The file {transfer_list_file} exists, to force the transfer delete it or use --force option")
                        if not args.dry_run:
                            with open(transfer_list_file, "a", encoding="utf-8") as transfer_files:
                                transfer_files.write(f"{abacus_bam} " + os.path.join(beluga_raw_reads_dir, f"{sample}.bam") + "\n")
                        else:
                            print(f"{abacus_bam} " + os.path.join(beluga_raw_reads_dir, f"{sample}.bam") + "\n")
                        update_db(connection,
                            sample,
                            readset,
                            library_type,
                            run_type,
                            run,
                            lane,
                            adapter1,
                            adapter2,
                            quality_offset,
                            bed,
                            fastq1,
                            fastq2,
                            readset_bam,
                            raw_reads_count,
                            raw_mean_coverage,
                            raw_median_insert_size,
                            raw_mean_insert_size,
                            raw_duplication_rate,
                            sample[-2:],
                            run_name)

                elif sample.endswith('RT'):
                    project_id = line_splitted[4]
                    # print(glob.glob(os.path.join(os.path.dirname(run_csv), f"Unaligned.{lane}", f"Project_{project_id}", f"Sample_{sample}_{library_id}", f"{sample}_{library_id}_S*_L00{lane}_R1_001.fastq.gz")))
                    # break
                    abacus_fastq1_file = os.path.join(os.path.dirname(args.run_csv), f"Unaligned.{lane}", f"Project_{project_id}", f"Sample_{sample}_{library_id}", f"{sample}_{library_id}_S*_L00{lane}_R1_001.fastq.gz")
                    abacus_fastq1 = "".join(glob.glob(abacus_fastq1_file))
                    abacus_fastq2_file = os.path.join(os.path.dirname(args.run_csv), f"Unaligned.{lane}", f"Project_{project_id}", f"Sample_{sample}_{library_id}", f"{sample}_{library_id}_S*_L00{lane}_R2_001.fastq.gz")
                    abacus_fastq2 = "".join(glob.glob(abacus_fastq2_file))
                    readset_bam = None
                    fastq1 = os.path.join("raw_reads", sample, f"{sample}_R1.fastq.gz")
                    fastq2 = os.path.join("raw_reads", sample, f"{sample}_R2.fastq.gz")
                    if not os.path.isfile(abacus_fastq1):
                        print(f"Warning: fastq file {abacus_fastq1_file} doesn't exist for sample {sample}.")
                        fastq1 = None
                    if not os.path.isfile(abacus_fastq2):
                        print(f"Warning: fastq file {abacus_fastq2_file} doesn't exist for sample {sample}.")
                        fastq2 = None
                    if fastq1 and fastq2:
                        # if os.path.getsize(transfer_list_file) != 0 and not args.force:
                        #     raise FileExistsError(f"The file {transfer_list_file} exists, to force the transfer delete it or use --force option")
                        if not args.dry_run:
                            with open(transfer_list_file, "a", encoding="utf-8") as transfer_files:
                                transfer_files.write(f"{abacus_fastq1} " + os.path.join(beluga_raw_reads_dir, f"{sample}_R1.fastq.gz") + "\n")
                                transfer_files.write(f"{abacus_fastq2} " + os.path.join(beluga_raw_reads_dir, f"{sample}_R2.fastq.gz") + "\n")
                        else:
                            print(f"{abacus_fastq1} " + os.path.join(beluga_raw_reads_dir, f"{sample}_R1.fastq.gz") + "\n")
                            print(f"{abacus_fastq2} " + os.path.join(beluga_raw_reads_dir, f"{sample}_R2.fastq.gz") + "\n")
                        update_db(connection,
                            sample,
                            readset,
                            library_type,
                            run_type,
                            run,
                            lane,
                            adapter1,
                            adapter2,
                            quality_offset,
                            bed,
                            fastq1,
                            fastq2,
                            readset_bam,
                            raw_reads_count,
                            raw_mean_coverage,
                            raw_median_insert_size,
                            raw_mean_insert_size,
                            raw_duplication_rate,
                            sample[-2:],
                            run_name)
                progress.update(index)

    connection.commit()
    connection.close()

    # Rsyncing back database to beluga:
    subprocess.run(f"""rsync -u {abacus_db} {args.user}@beluga.calculcanada.ca:{beluga_db}""", check=True, shell=True)

    # Doing Globus transfer
    if not args.dry_run:
        print("Submitting globus transfer...")
        subprocess.run(f"""
    module load mugqic/globus-cli/3.7.0
    sub_id="$(globus task generate-submission-id)"
    globus transfer -verify-checksum --submission-id $sub_id --label "{run_name}" --batch {transfer_list_file} {abacus_endpoint} {beluga_endpoint}""",
            check=True,
            shell=True
            )

def update_db(connection,
    sample,
    readset,
    library_type,
    run_type,
    run,
    lane,
    adapter1,
    adapter2,
    quality_offset,
    bed,
    fastq1,
    fastq2,
    readset_bam,
    raw_reads_count,
    raw_mean_coverage,
    raw_median_insert_size,
    raw_mean_insert_size,
    raw_duplication_rate,
    sample_type,
    run_name):
    # Updating readset table
    update_readset_table(
        connection,
        sample,
        readset,
        library_type,
        run_type,
        run,
        lane,
        adapter1,
        adapter2,
        quality_offset,
        bed,
        fastq1,
        fastq2,
        readset_bam
        )
    # Updating key metric table
    fails = []
    flags = []
    raw_reads_count_check(raw_reads_count, sample_type, fails, flags)
    raw_mean_coverage_check(raw_mean_coverage, sample_type, fails)
    raw_duplication_rate_check(raw_duplication_rate, sample_type, fails, flags)
    try:
        flags.extend(extract_value(connection, "key_metric", sample, "flag").split(";"))
        fails.extend(extract_value(connection, "key_metric", sample, "fail").split(";"))
    except AttributeError:
        pass
    flags = ';'.join(set(flags))
    fails = ';'.join(set(fails))
    update_key_metric_table_run_metrics(
        connection,
        sample,
        raw_reads_count,
        raw_mean_coverage,
        raw_median_insert_size,
        raw_mean_insert_size,
        raw_duplication_rate,
        flags,
        fails
    )
    # Updating sample table
    result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
    patient = result.group(1)
    cohort = result.group(2)
    institution = result.group(3)
    update_sample_table(connection, sample, patient, run_name)
    # Updating file_location table
    if sample_type == "DT":
        dna_abacus_raw_bam_t = readset_bam
        dna_abacus_raw_bam_n = None
    elif sample_type == "DN":
        dna_abacus_raw_bam_n = readset_bam
        dna_abacus_raw_bam_t = None
    else:
        dna_abacus_raw_bam_n = None
        dna_abacus_raw_bam_t = None
    update_fileloc_details_run_processing(
        patient,
        dna_abacus_raw_bam_t=dna_abacus_raw_bam_t,
        dna_abacus_raw_bam_n=dna_abacus_raw_bam_n,
        rna_abacus_raw_fastq1=fastq1,
        rna_abacus_raw_fastq2=fastq2
    )

if __name__ == '__main__':
    main()
