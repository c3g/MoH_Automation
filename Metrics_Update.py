#!/usr/bin/env python3

import glob
import os
import re
import csv
import h5py
import numpy as np
import sys
import progressbar
from  DB_OPS import update_metrics_db,create_connection,extract_sample_field,extract_sample_names, extract_fileloc_field, extract_value

WIDGETS = [' [', progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') - ', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']


def main():
    connection = create_connection("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db")
    # connection = create_connection("/scratch/stretenp/moh_test/MOH_analysis.db")

    patients = extract_sample_names(connection)
    #TEST CASE
    #Samples = ['MoHQ-JG-9-10','MoHQ-MU-12-122','MoHQ-MU-12-1224']

    samples_list = []
    paired_samples_dict = {}

    # All_Samples = []
    # PAIRED_SAMPLES = dict()
    print("Fetching Database by Patient...")
    with progressbar.ProgressBar(max_value=len(patients), widgets=WIDGETS) as progress:
        for index, patient in enumerate(patients, 1):
            # run_sample = extract_sample_field(connection, patient, "run")
            dna_n_sample = extract_sample_field(connection, patient, "DNA_N")
            if dna_n_sample not in ("NA", ""):
                samples_list.append(dna_n_sample)
                try:
                    run_sample = extract_fileloc_field(connection, patient, "Run_Proc_BAM_DNA_N").split("/")[7]
                except IndexError:
                    run_sample = ""
                paired_samples_dict[dna_n_sample] = (patient, run_sample)
            dna_t_sample = extract_sample_field(connection, patient, "DNA_T")
            if dna_t_sample not in ("NA", ""):
                samples_list.append(dna_t_sample)
                try:
                    run_sample = extract_fileloc_field(connection, patient, "Run_Proc_BAM_DNA_T").split("/")[7]
                except IndexError:
                    run_sample = ""
                paired_samples_dict[dna_t_sample] = (patient, run_sample)
            rna_sample = extract_sample_field(connection, patient, "RNA")
            if rna_sample not in ("NA", ""):
                samples_list.append(rna_sample)
                try:
                    run_sample = extract_fileloc_field(connection, patient, "Run_Proc_fastq_1_RNA").split("/")[7]
                except IndexError:
                    run_sample = ""
                paired_samples_dict[rna_sample] = (patient, run_sample)
            progress.update(index)
    extract_data(samples_list, connection, paired_samples_dict)
    print("Committing changes to Database...")
    connection.commit()
    connection.close()
    print("...Done.")

def extract_data(samples_list, connection, paired_samples_dict):
    print("Updating metrics in Database by Sample...")
    with progressbar.ProgressBar(max_value=len(samples_list), widgets=WIDGETS) as progress:
        for index, sample in enumerate(samples_list, 1):
            patient = paired_samples_dict[sample][0]
            run = paired_samples_dict[sample][1]
            if sample.endswith('DN'):
                sample_type = 'DN'
            elif sample.endswith('DT'):
                sample_type = 'DT'
            elif sample.endswith('RT'):
                sample_type = 'RT'
            flags = []
            fails = []

            # WGS_Bases_Over_Q30
            dna_bases_over_q30_percent = extract_bs_over_q30(sample, sample_type)
            try:
                if int(dna_bases_over_q30_percent)<75 and sample_type in ('DN', 'DT'):
                    fails.append('WGS_Bases_Over_Q30')
                elif int(dna_bases_over_q30_percent)<80 and sample_type in ('DN', 'DT'):
                    flags.append('WGS_Bases_Over_Q30')
            except (TypeError, ValueError):
                pass

            # WGS_Min_Aligned_Reads_Delivered
            dna_aligned_reads_count = extract_min_aln_rds(sample, patient)
            # QC Gates
            # try:
            #     if int(dna_aligned_reads_count)<260000000 and sample_type == 'DN':
            #         fails.append('WGS_Min_Aligned_Reads_Delivered')
            #     elif int(dna_aligned_reads_count)<660000000 and sample_type == 'DN':
            #         flags.append('WGS_Min_Aligned_Reads_Delivered')
            #     elif int(dna_aligned_reads_count)<530000000 and sample_type == 'DT':
            #         fails.append('WGS_Min_Aligned_Reads_Delivered')
            #     elif int(dna_aligned_reads_count)<1330000000 and sample_type == 'DT':
            #         flags.append('WGS_Min_Aligned_Reads_Delivered')
            # except (TypeError, ValueError):
            #     pass

            # Raw_Reads_Count, Raw_Mean_Coverage, Raw_Median_Insert_Size, Raw_Mean_Insert_Size, Raw_Duplication_Rate
            raw_reads_count, raw_mean_coverage, raw_median_insert_size, raw_mean_insert_size, raw_duplication_rate = parse_run_metrics(sample, run)
            # QC Gates
            # try:
            #     if float(raw_mean_coverage)<30 and sample_type == 'DN':
            #         fails.append('Raw_Mean_Coverage')
            #     elif float(raw_mean_coverage)<80 and sample_type == 'DT':
            #         fails.append('Raw_Mean_Coverage')
            # except (TypeError, ValueError):
            #     pass
            try:
                if float(raw_reads_count)<80000000 and sample_type == 'RT':
                    fails.append('Raw_Reads_Count')
                elif float(raw_reads_count)<100000000 and sample_type == 'RT':
                    flags.append('Raw_Reads_Count')
            except (TypeError, ValueError):
                pass
            # QC Gates
            # try:
            #     if float(raw_duplication_rate)>50 and sample_type in ('DT', 'DN'):
            #         fails.append('Raw_Duplication_Rate')
            #     elif float(raw_duplication_rate)>20 and sample_type in ('DT', 'DN'):
            #         flags.append('Raw_Duplication_Rate')
            # except (TypeError, ValueError):
            #     pass

            # WGS_Dedup_Coverage
            dna_dedup_coverage = extract_dedup_coverage(sample)
            try:
                if float(dna_dedup_coverage)<30 and sample_type == 'DN':
                    fails.append('WGS_dedup_coverage')
                elif float(dna_dedup_coverage)<80 and sample_type == 'DT':
                    fails.append('WGS_dedup_coverage')
            except (TypeError, ValueError):
                pass

            # Median_Insert_Size, Mean_Insert_Size
            median_insert_size, mean_insert_size = extract_insert_size(sample, patient, sample_type)
            try:
                if float(median_insert_size)<300:
                    flags.append('Median_Insert_Size')
                elif float(median_insert_size)<150:
                    fails.append('Median_Insert_Size')
            except (TypeError, ValueError):
                pass

            # WGS_Contamination
            dna_contamination = extract_contamination(patient, sample_type)
            try:
                if float(dna_contamination)>0.5:
                    fails.append('WGS_Contamination')
                elif float(dna_contamination)>0.05:
                    flags.append('WGS_Contamination')
            except (TypeError, ValueError):
                pass

            # Concordance
            dna_concordance = extract_concordance(patient, sample, sample_type)
            try:
                if float(dna_concordance)<99:
                    fails.append('Concordance')
            except (TypeError, ValueError):
                pass

            # Purity
            dna_tumour_purity = extract_purity(sample, patient)
            try:
                if float(dna_tumour_purity)<30:
                    fails.append('Purity')
            except (TypeError, ValueError):
                pass

            # WTS_Aligned_Reads
            rna_aligned_reads_count = extract_rna_aligned_reads_count(sample)


            # WTS_Expression_Profiling_Efficiency, WTS_rRNA_contamination
            rna_expression_profiling_efficiency, rna_ribosomal_contamination_count = parse_rnaseqc2(sample)
            try:
                if float(rna_ribosomal_contamination_count)>0.35:
                    fails.append('WTS_rRNA_contamination')
                elif float(rna_ribosomal_contamination_count)>0.1:
                    flags.append('WTS_rRNA_contamination')
            except (TypeError, ValueError):
                pass
            # # WTS_Expression_Profiling_Efficiency
            # rna_expression_profiling_efficiency = extract_rna_expression_profiling_efficiency(sample)

            # # WTS_rRNA_contamination
            # rrna_count = extract_rna_ribosomal(sample)
            # try:
            #     rna_ribosomal_contamination_count = int(rrna_count)/int(rna_aligned_reads_count)
            #     if float(rna_ribosomal_contamination_count)>0.35:
            #         fails.append('WTS_rRNA_contamination')
            #     elif float(rna_ribosomal_contamination_count)>0.1:
            #         flags.append('WTS_rRNA_contamination')
            # except (TypeError, ValueError):
            #     rna_ribosomal_contamination_count = "NA"

            # Flags
            tmp_flags = extract_value(connection, "KEY_METRICS", sample, "Flags").split(";")
            if len(tmp_flags) == 1 and tmp_flags[0] == "NA":
                pass
            else:
                flags.extend(tmp_flags)
            tmp_fails = extract_value(connection, "KEY_METRICS", sample, "Fails").split(";")
            if len(tmp_fails) == 1 and tmp_fails[0] == "NA":
                pass
            else:
                fails.extend(tmp_fails)

            try:
                flags.remove("NA")
            except ValueError:
                pass
            flags = ';'.join(set(flags))
            if not flags:
                flags = "NA"
            try:
                fails.remove("NA")
            except ValueError:
                pass
            fails = ';'.join(set(fails))
            if not fails:
                fails = "NA"

            update_metrics_db(
                conn=connection,
                sample=sample,
                WGS_Bases_Over_Q30=dna_bases_over_q30_percent,
                WGS_Min_Aligned_Reads_Delivered=dna_aligned_reads_count,
                Raw_Mean_Coverage=raw_mean_coverage,
                WGS_Dedup_Coverage=dna_dedup_coverage,
                Median_Insert_Size=median_insert_size,
                Mean_Insert_Size=mean_insert_size,
                Raw_Duplication_Rate=raw_duplication_rate,
                WGS_Contamination=dna_contamination,
                Raw_Reads_Count=raw_reads_count,
                WTS_Aligned_Reads=rna_aligned_reads_count,
                WTS_rRNA_contamination=rna_ribosomal_contamination_count,
                WTS_Expression_Profiling_Efficiency=rna_expression_profiling_efficiency,
                Concordance=dna_concordance,
                Purity=dna_tumour_purity,
                Flags=flags,
                Fails=fails,
                Raw_Median_Insert_Size=raw_median_insert_size,
                Raw_Mean_Insert_Size=raw_mean_insert_size
                )
            progress.update(index)


def parse_run_metrics(sample, run):
    """Raw_Reads_Count, Raw_Mean_Coverage, Raw_Median_Insert_Size, Raw_Mean_Insert_Size, Raw_Duplication_Rate"""
    raw_reads_count = "NA"
    raw_mean_coverage = "NA"
    raw_median_insert_size = "NA"
    raw_mean_insert_size = "NA"
    raw_duplication_rate = "NA"
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics', run + '-run.align_bwa_mem.csv')
        with open(filename, 'r', encoding="utf-8") as file:
            for line in file:
                parsed_line = line.split(",")
                if sample == parsed_line[6]:
                    if parsed_line[12]:
                        raw_reads_count = parsed_line[12]
                    if parsed_line[41]:
                        raw_mean_coverage = parsed_line[41]
                    if parsed_line[37]:
                        raw_median_insert_size = parsed_line[37]
                    if parsed_line[38]:
                        raw_mean_insert_size = parsed_line[38]
                    if parsed_line[15]:
                        raw_duplication_rate = parsed_line[15]
    except FileNotFoundError:
        raw_reads_count = "NA"
        raw_mean_coverage = "NA"
        raw_median_insert_size = "NA"
        raw_mean_insert_size = "NA"
        raw_duplication_rate = "NA"
    return raw_reads_count, raw_mean_coverage, raw_median_insert_size, raw_mean_insert_size, raw_duplication_rate


def parse_rnaseqc2(sample):
    """WTS_Expression_Profiling_Efficiency, WTS_rRNA_contamination"""
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample, 'rnaseqc2', sample + '.sorted.mdup.bam.metrics.tsv')
        with open(filename, 'r', encoding="utf-8") as file:
            for line in file:
                if line.startswith("Expression Profiling Efficiency"):
                    rna_expression_profiling_efficiency = line.split("\t")[1]
                elif line.startswith("rRNA Rate"):
                    rna_ribosomal_contamination_count = line.split("\t")[1]
    except FileNotFoundError:
        rna_expression_profiling_efficiency = rna_ribosomal_contamination_count = "NA"
    return rna_expression_profiling_efficiency, rna_ribosomal_contamination_count

def extract_rna_ribosomal(sample):
    """WTS_rRNA_contamination"""
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample, 'rnaseqc', sample, sample + '.rRNA_counts.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            line = lines[0]
            fields = line.split("\t")
            ret = fields[0]
    except FileNotFoundError:
        ret = "NA"
    return ret


def extract_rna_expression_profiling_efficiency(sample):
    """WTS_Expression_Profiling_Efficiency"""
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample, 'rnaseqc2', sample + '.sorted.mdup.bam.metrics.tsv')
        with open(filename, 'r', encoding="utf-8") as file:
            for line in file:
                if line.startswith("Expression Profiling Efficiency"):
                    rna_expression_profiling_efficiency = line.split("\t")[1]
    except FileNotFoundError:
        rna_expression_profiling_efficiency = "NA"
    return rna_expression_profiling_efficiency


def extract_rna_aligned_reads_count(sample):
    """WTS_Aligned_Reads"""
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics', sample, sample + '.alignment_summary_metrics')
        with open(filename, 'r', encoding="utf-8") as file:
            for line in file:
                if line.startswith("PAIR"):
                    rna_aligned_reads_count = line.split("\t")[6]
    except FileNotFoundError:
        rna_aligned_reads_count = "NA"
    return rna_aligned_reads_count


def extract_purity(sample, patient):
    """Purity"""
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, 'purple', sample + '.purple.purity.tsv')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            line = lines[1]
            fields = line.split("\t")
            ret = float(fields[0])*100
    except FileNotFoundError:
        ret = "NA"
    return ret


def extract_contamination(patient, sample_type):
    """WGS_Contamination"""
    ret = "NA"
    if sample_type in ('DN', 'DT'):
        # The file is named after tumour sample only
        filename = glob.glob(os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics', patient + '-*DT.contamination.tsv'))
        # Test if unsure about finding more than 1 concordance file for normal sample based on the glob above.
        # It has to print nothing to be ok, otherwise it means a manual check is required.
        if len(filename) > 1:
            print(f" WARNING: Manual check required for patient {patient} as more than 1 concordance file is found (by default the first one is selected for tghe metric): {filename}")
        try:
            with open(filename[0], 'r', encoding="utf-8") as file:
                for line in file:
                    if line.startswith('Normal') and sample_type == 'DN':
                        ret = line.split(" ")[-1][:-2]
                    elif line.startswith('Tumor') and sample_type == 'DT':
                        ret = line.split(" ")[-1][:-2]
        except (FileNotFoundError, IndexError):
            ret = "NA"
    return ret


def extract_concordance(patient, sample, sample_type):
    """Concordance"""
    ret = "NA"
    if sample_type in ('DN', 'DT'):
        # The file is named after tumour sample only
        filename = glob.glob(os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics', f"{patient}-*DT.concordance.tsv"))
        # Test if unsure about finding more than 1 concordance file for normal sample based on the glob above.
        # It has to print nothing to be ok, otherwise it means a manual check is required.
        if len(filename) > 1:
            print(f" WARNING: Manual check required for patient {patient} as more than 1 concordance file is found(by default the first one is selected for tghe metric): {filename}")
        try:
            with open(filename[0], 'r', encoding="utf-8") as file:
                for line in file:
                    if line.startswith('Concordance'):
                        ret = line.split(" ")[-1][:-2]
        except (FileNotFoundError, IndexError):
            ret = "NA"
    return ret


def extract_insert_size(sample, patient, sample_type):
    """Median_Insert_Size, Mean_Insert_Size"""
    median_insert_size = mean_insert_size = "NA"
    try:
        if sample_type == 'RT':
            filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/kallisto', sample, 'abundance.h5')
            with h5py.File(filename, 'r') as h5file:
                try:
                    frequencies = np.asarray(h5file['aux']['fld'])
                except OSError:
                    raise Exception(f"ERROR: file {filename} for sample {sample} is corrupted.")
                values = np.arange(0, len(frequencies), 1)
                # Calculus coming from https://stackoverflow.com/questions/46086663/how-to-get-mean-and-standard-deviation-from-a-frequency-distribution-table
                ord = np.argsort(values)
                cdf = np.cumsum(frequencies[ord])
                median_insert_size = values[ord][np.searchsorted(cdf, cdf[-1] // 2)].astype('float64')
                mean_insert_size = np.around(np.average(values, weights=frequencies), decimals=1)
        elif sample_type in ('DN', 'DT'):
            filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'qualimap', sample, "raw_data_qualimapReport/insert_size_histogram.txt")
            dataset = np.genfromtxt(fname=filename, delimiter="\t", skip_header=1)
            frequencies = dataset[:, 1]
            values = dataset[:, 0]
            # Calculus coming from https://stackoverflow.com/questions/46086663/how-to-get-mean-and-standard-deviation-from-a-frequency-distribution-table
            ord = np.argsort(values)
            cdf = np.cumsum(frequencies[ord])
            median_insert_size = values[ord][np.searchsorted(cdf, cdf[-1] // 2)]
            mean_insert_size = np.around(np.average(values, weights=frequencies), decimals=1)
    except FileNotFoundError:
        pass
    return median_insert_size, mean_insert_size


def extract_dedup_coverage(sample):
    """WGS_Dedup_Coverage"""
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'qualimap', sample, 'genome_results.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            line = lines[71]
            # line is mean coverageData = 304.9902X
            metrics = line.split(" ")
            ret = float(metrics[-1].replace('X', ''))
    except (FileNotFoundError, ValueError):
        ret = "NA"
    return ret


def extract_min_aln_rds(sample, patient):
    """WGS_Min_Aligned_Reads_Delivered"""
    ret = "NA"
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', patient + '.multiqc_data', 'multiqc_general_stats.txt')
        with open(filename, 'r', encoding="utf-8") as csvfile:
            reader = csv.DictReader(csvfile, delimiter="\t")
            for row in reader:
                if row["Sample"] == sample and row['QualiMap_mqc-generalstats-qualimap-mapped_reads']:
                    ret = row["QualiMap_mqc-generalstats-qualimap-mapped_reads"]
    except (FileNotFoundError, KeyError):
        ret = "NA"
    return ret


def extract_bs_over_q30(sample, sample_type):
    """WGS_Bases_Over_Q30"""
    try:
        if sample_type in ('DT', 'DN'):
            filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'picard_metrics', sample + '.all.metrics.quality_distribution_metrics')
        elif sample_type == 'RT':
            filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample + '.quality_distribution_metrics')
        tester = re.compile('(\d+)\W+(\d+)')
        with open(filename, 'r', encoding="utf-8") as file:
            above_30 = 0
            below_30 = 0
            for line in file:
                if line[:1].isdigit():
                    test = tester.match(line)
                    qual = test.group(1)
                    count = test.group(2)
                    if int(qual) < 30:
                        below_30 += int(count)
                    else :
                        above_30 += int(count)
            percent_abv = round((above_30/(above_30+below_30))*100, 2)
    except FileNotFoundError:
        percent_abv = "NA"
    return percent_abv


if __name__ == '__main__':
    main()
    #Update db with the objects
