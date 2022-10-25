#!/usr/bin/env python3

import glob
import os
import re
import progressbar
from  DB_OPS import update_metrics_db,create_connection,extract_sample_field,extract_sample_names, extract_fileloc_field, extract_value
#TO DO::
#FIX!

WIDGETS = [' [', progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') - ', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']


def main():
    # widgets = [' [', progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') - ', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']
    connection = create_connection(r"/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db")

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
            if dna_n_sample:
                samples_list.append(dna_n_sample)
                try:
                    run_sample = extract_fileloc_field(connection, patient, "Run_Proc_BAM_DNA_N").split("/")[7]
                except IndexError:
                    run_sample = ""
                paired_samples_dict[dna_n_sample] = (patient, run_sample)
            dna_t_sample = extract_sample_field(connection, patient, "DNA_T")
            if dna_t_sample:
                samples_list.append(dna_t_sample)
                try:
                    run_sample = extract_fileloc_field(connection, patient, "Run_Proc_BAM_DNA_T").split("/")[7]
                except IndexError:
                    run_sample = ""
                paired_samples_dict[dna_t_sample] = (patient, run_sample)
            rna_sample = extract_sample_field(connection, patient, "RNA")
            if rna_sample:
                samples_list.append(rna_sample)
                try:
                    run_sample = extract_fileloc_field(connection, patient, "Run_Proc_fastq_1_RNA").split("/")[7]
                except IndexError:
                    run_sample = ""
                paired_samples_dict[rna_sample] = (patient, run_sample)
            progress.update(index)
    # with progressbar.ProgressBar(max_value=len(patients), widgets=widgets) as bar:
    #     for index, sample in enumerate(patients, 1):
    #         dna_n = extract_sample_field(connection,sample,"DNA_N")
    #         if dna_n != "NA" and dna_n != "":
    #             All_Samples.append(dna_n)
    #             PAIRED_SAMPLES[dna_n] = sample
    #         dna_t = extract_sample_field(connection,sample,"DNA_T")
    #         if dna_t != "NA" and dna_t != "":
    #             All_Samples.append(dna_t)
    #             PAIRED_SAMPLES[dna_t] = sample
    #         rna = extract_sample_field(connection,sample,"RNA")
    #         if rna != "NA" and rna != "":
    #             All_Samples.append(rna)
    #             PAIRED_SAMPLES[rna] = sample
    #         bar.update(index)
    #print(All_Samples)
    #All_Samples = [x for x in All_Samples if x != "NA" and x != ""]
    #print(len(All_Samples))
    #print(len(PAIRED_SAMPLES))
    #exit()
    extract_data(samples_list, connection, paired_samples_dict)
    # extract_data(All_Samples, connection, PAIRED_SAMPLES)
    print("Committing changes to Database...")
    connection.commit()
    connection.close()
    print("Done.")

def extract_data(samples_list, connection, paired_samples_dict):
    print("Updating metrics in Database by Sample...")
    # with progressbar.ProgressBar(max_value=len(SAMP), widgets=WIDGETS) as bar:
    #     for index, Sample in enumerate(SAMP, 1):
    #         PATIENT = PAIRED_SAMPLES[Sample]
    #         TYPE=''
    #         if Sample.endswith('N'):
    #             TYPE = 'N'
    #         elif Sample.endswith('T'):
    #             TYPE = 'T'
    #         flag=[]
    #         fail=[]
    with progressbar.ProgressBar(max_value=len(samples_list), widgets=WIDGETS) as progress:
        for index, sample in enumerate(samples_list, 1):
            # print(sample + "\n")
            patient = paired_samples_dict[sample][0]
            run = paired_samples_dict[sample][1]
            # sample_type = ''
            if sample.endswith('DN'):
                sample_type = 'DN'
            elif sample.endswith('DT'):
                sample_type = 'DT'
            elif sample.endswith('RT'):
                sample_type = 'RT'
            flags = []
            fails = []

            dna_bases_over_q30_percent = extract_bs_over_q30(sample, sample_type)
            try:
                if int(dna_bases_over_q30_percent)<75 and sample_type in ('DN', 'DT'):
                    fails.append('WGS_Bases_Over_Q30')
                elif int(dna_bases_over_q30_percent)<80 and sample_type in ('DN', 'DT'):
                    flags.append('WGS_Bases_Over_Q30')
            except TypeError:
                pass
            # if dna_bases_over_q30_percent is None:
            #     dna_bases_over_q30_percent = ''
            # elif float(dna_bases_over_q30_percent)<75:
            #     fail.append('WGS_Bases_Over_Q30')
            # elif float(dna_bases_over_q30_percent)<80:
            #     flag.append('WGS_Bases_Over_Q30')

            dna_aligned_reads_count = extract_min_aln_rds(sample, patient)
            try:
                if int(dna_aligned_reads_count)<260000000 and sample_type == 'DN':
                    fails.append('WGS_Min_Aligned_Reads_Delivered')
                elif int(dna_aligned_reads_count)<660000000 and sample_type == 'DN':
                    flags.append('WGS_Min_Aligned_Reads_Delivered')
                elif int(dna_aligned_reads_count)<530000000 and sample_type == 'DT':
                    fails.append('WGS_Min_Aligned_Reads_Delivered')
                elif int(dna_aligned_reads_count)<1330000000 and sample_type == 'DT':
                    flags.append('WGS_Min_Aligned_Reads_Delivered')
            except TypeError:
                pass
            # if dna_aligned_reads_count in (None, ''):
            #     dna_aligned_reads_count = ''
            # elif (float(dna_aligned_reads_count)<260000000 and sample_type == 'DN'):
            #     fail.append('WGS_Min_Aligned_Reads_Delivered')
            # elif (float(dna_aligned_reads_count)<660000000 and sample_type == 'DN'):
            #     flag.append('WGS_Min_Aligned_Reads_Delivered')
            # elif (float(dna_aligned_reads_count)<530000000 and sample_type == 'DT'):
            #     fail.append('WGS_Min_Aligned_Reads_Delivered')
            # elif (float(dna_aligned_reads_count)<1330000000 and sample_type == 'DT'):
            #     flag.append('WGS_Min_Aligned_Reads_Delivered')
            # if WGS_Min_Aligned_Reads_Delivered == None or WGS_Min_Aligned_Reads_Delivered == '' :
            #     WGS_Min_Aligned_Reads_Delivered = ''
            # elif (float(WGS_Min_Aligned_Reads_Delivered)<260000000 and TYPE == 'N'):
            #     fail.append('WGS_Min_Aligned_Reads_Delivered')
            # elif (float(WGS_Min_Aligned_Reads_Delivered)<660000000 and TYPE == 'N'):
            #     flag.append('WGS_Min_Aligned_Reads_Delivered')
            # elif (float(WGS_Min_Aligned_Reads_Delivered)<530000000 and TYPE == 'T'):
            #     fail.append('WGS_Min_Aligned_Reads_Delivered')
            # elif (float(WGS_Min_Aligned_Reads_Delivered)<1330000000 and TYPE == 'T'):
            #     flag.append('WGS_Min_Aligned_Reads_Delivered')


            raw_reads_count, raw_mean_coverage, raw_median_insert_size, raw_mean_insert_size, raw_duplication_rate = parse_run_metrics(sample, run)
            try:
                if float(raw_mean_coverage)<30 and sample_type == 'DN':
                    fails.append('Raw_Mean_Coverage')
                elif float(raw_mean_coverage)<80 and sample_type == 'DT':
                    fails.append('Raw_Mean_Coverage')
            except TypeError:
                pass
            try:
                if float(raw_reads_count)<80000000 and sample_type == 'RT':
                    fails.append('Raw_Reads_Count')
                elif float(raw_reads_count)<100000000 and sample_type == 'RT':
                    flags.append('Raw_Reads_Count')
            except TypeError:
                pass
            try:
                if float(raw_duplication_rate)>50 and sample_type in ('DT', 'DN'):
                    fails.append('Raw_Duplication_Rate')
                elif float(raw_duplication_rate)>20 and sample_type in ('DT', 'DN'):
                    flags.append('Raw_Duplication_Rate')
            except TypeError:
                pass

            # dna_raw_coverage = extract_raw_coverage(sample)
            # if  dna_raw_coverage in (None, ''):
            #     dna_raw_coverage = ''
            # elif (float(dna_raw_coverage)<30 and sample_type == 'N'):
            #     fail.append('WGS_Raw_Coverage')
            # elif (float(dna_raw_coverage)<80 and sample_type == 'T'):
            #     fail.append('WGS_Raw_Coverage')

            # #WGS_Duplication_Rate
            # WGS_duplicates = extract_sambama_dups(Sample,PATIENT)
            # if  WGS_duplicates == None or WGS_duplicates == '' :
            #     WGS_duplicates = ''
            # elif (float(WGS_duplicates)>50):
            #     fail.append('WGS_Duplication_Rate')
            # elif (float(WGS_duplicates)>20):
            #     flag.append('WGS_Duplication_Rate')

            # WGS_Raw_Coverage = extract_raw_coverage(Sample)
            # if  WGS_Raw_Coverage == None or WGS_Raw_Coverage == '' :
            #     WGS_Raw_Coverage = ''
            # elif (float(WGS_Raw_Coverage)<30 and TYPE == 'N'):
            #     fail.append('WGS_Raw_Coverage')
            # elif (float(WGS_Raw_Coverage)<80 and TYPE == 'T'):
            #     fail.append('WGS_Raw_Coverage')


            dna_dedup_coverage = extract_dedup_coverage(sample)
            # WGS_Dedup_Coverage = extract_ded_coverage(Sample,PATIENT)
            # if WGS_Dedup_Coverage  == None or WGS_Dedup_Coverage == '' :
            #     WGS_Dedup_Coverage = ''

            median_insert_size = extract_insert_size(sample, patient, sample_type)
            try:
                if float(median_insert_size)<300:
                    flags.append('Median_Insert_Size')
                elif float(median_insert_size)<150:
                    fails.append('Median_Insert_Size')
            except TypeError:
                pass
            # Median_Insert_Size = extract_insert_size(Sample,PATIENT)
            # if  Median_Insert_Size == None or Median_Insert_Size == '' :
            #     Median_Insert_Size = ''
            # elif (float(Median_Insert_Size)<300):
            #     flag.append('Median_Insert_Size')
            # elif (float(Median_Insert_Size)<150):
            #     fail.append('Median_Insert_Size')

            #WGS_Contamination
            dna_contamination = extract_contamination(patient, sample_type)
            try:
                if float(dna_contamination)>5:
                    fails.append('WGS_Contamination')
            except TypeError:
                pass
            # WGS_Contamination = extract_contamination(Sample,PATIENT)
            # if WGS_Contamination  == None or WGS_Contamination == '' :
            #     WGS_Contamination = ''
            # elif (float(WGS_Contamination)>5):
            #     fail.append('WGS_Contamination')

            #Concordance
            dna_concordance = extract_concordance(patient, sample_type)
            try:
                if float(dna_concordance)<99:
                    fails.append('Concordance')
            except TypeError:
                pass
            # Concordance = extract_concordance(Sample,PATIENT)
            # if  Concordance == None or Concordance == '' :
            #     Concordance = ''
            # elif (float(Concordance)<99):
            #     fail.append('Concordance')

            # Tumor_Purity
            dna_tumour_purity = extract_purity(sample, patient)
            try:
                if float(dna_tumour_purity)<30:
                    fails.append('Purity')
            except TypeError:
                pass
            # Purity = extract_purity(Sample,PATIENT)
            # if  Purity == None or Purity == '' :
            #     Purity = ''
            # elif (float(Purity)<30):
            #     fail.append('Purity')

            #WTS_Clusters
            # WTS_Clusters = extract_WTS_Clusters(Sample)
            # if  WTS_Clusters == None or WTS_Clusters == '' :
            #    WTS_Clusters  = ''
            # elif (float(WTS_Clusters)<80000000):
            #     fail.append('WTS_Clusters')
            # elif (float(WTS_Clusters)<100000000):
            #     flag.append('WTS_Clusters')

            #WTS_Exonic_Rate
            rna_aligned_reads_count, rna_exonic_rate = parse_rnaseqc_metrics_tmp(sample)
            try:
                if float(rna_exonic_rate)<0.6:
                    fails.append('WTS_Exonic_Rate')
                elif float(rna_exonic_rate)<0.8:
                    flags.append('WTS_Exonic_Rate')
            except TypeError:
                pass
            # WTS_Exonic_Rate = extract_WTS_exonic(Sample)
            # if WTS_Exonic_Rate  == None or WTS_Exonic_Rate == '' :
            #     WTS_Exonic_Rate = ''
            # elif (float(WTS_Exonic_Rate)<0.6):
            #     fail.append('WTS_Exonic_Rate')
            # elif (float(WTS_Exonic_Rate)<0.8):
            #     flag.append('WTS_Exonic_Rate')

            #WTS_Unique_Reads
            # WTS_Unique_Reads = extract_WTS_unique(Sample)
            # if  WTS_Unique_Reads == None or WTS_Unique_Reads == '' :
            #     WTS_Unique_Reads = ''

            # WTS_rRA_contamination
            rrna_count = extract_rna_ribosomal(sample)
            try:
                rna_ribosomal_contamination_count = int(rrna_count)/int(rna_aligned_reads_count)
                if float(rna_ribosomal_contamination_count)>0.35:
                    fails.append('WTS_rRNA_contamination')
                elif float(rna_ribosomal_contamination_count)>0.1:
                    flags.append('WTS_rRNA_contamination')
            except (TypeError, ValueError):
                rna_ribosomal_contamination_count = None
            # rRNA_count = extract_WTS_rRNA(Sample)
            # if  rRNA_count == None or rRNA_count == '' or WTS_Unique_Reads == '':
            #     WTS_rRNA_contamination = ''
            # else:
            #     WTS_rRNA_contamination = int(rRNA_count)/int(WTS_Unique_Reads)
            #     if (float(WTS_rRNA_contamination)>0.35):
            #         fail.append('WTS_rRNA_contamination')
            #     elif (float(WTS_rRNA_contamination)>0.1):
            #         flag.append('WTS_rRNA_contamination')

            # mean_ins_size = extract_mean_map_ins_size(Sample)
            # med_ins_size = extract_med_map_ins_size(Sample)

            # Flags
            flags.extend(extract_value(connection, "KEY_METRICS", sample, "Flags").split(";"))
            fails.extend(extract_value(connection, "KEY_METRICS", sample, "Fails").split(";"))
            flags = ';'.join(set(flags))
            fails = ';'.join(set(fails))
            # Yellow_Flags=';'.join(flag)
            # Red_Flags=';'.join(fail)

            #print(Sample,dna_bases_over_q30_percent,WGS_Min_Aligned_Reads_Delivered,WGS_Raw_Coverage,WGS_Dedup_Coverage,Median_Insert_Size,WGS_duplicates,WGS_Contamination,WTS_Clusters,WTS_Unique_Reads,WTS_Exonic_Rate,WTS_rRNA_contamination,Concordance,Purity,Yellow_Flags,Red_Flags,mean_ins_size,med_ins_size)
            # update_metrics_db(
            #     connection,
            #     Sample,
            #     WGS_Bases_Over_Q30,
            #     WGS_Min_Aligned_Reads_Delivered,
            #     WGS_Raw_Coverage,
            #     WGS_Dedup_Coverage,
            #     Median_Insert_Size,
            #     WGS_duplicates,
            #     WGS_Contamination,
            #     WTS_Clusters,
            #     WTS_Unique_Reads,
            #     WTS_Exonic_Rate,
            #     WTS_rRNA_contamination,
            #     Concordance,
            #     Purity,
            #     Yellow_Flags,
            #     Red_Flags,
            #     mean_ins_size,
            #     med_ins_size
            #     )

            update_metrics_db(
                connection,
                sample,
                dna_bases_over_q30_percent,
                dna_aligned_reads_count,
                raw_mean_coverage,
                dna_dedup_coverage,
                median_insert_size,
                raw_duplication_rate,
                dna_contamination,
                raw_reads_count,
                rna_aligned_reads_count,
                rna_exonic_rate,
                rna_ribosomal_contamination_count,
                dna_concordance,
                dna_tumour_purity,
                flags,
                fails,
                raw_mean_insert_size,
                raw_median_insert_size
                )
            progress.update(index)

            #print (Sample)
            #print (dna_bases_over_q30_percent)
            #print (WGS_Min_Aligned_Reads_Delivered)
            #print (WGS_Raw_Coverage)
            #print (WGS_Dedup_Coverage)
            #print (WGS_duplicates)
            #print (Median_Insert_Size)
            #print (mean_ins_size)
            #print (med_ins_size)
            #print (WGS_Contamination)
            #print (Concordance)
            #print (Purity)
            #print (WTS_Clusters)
            #print (WTS_Exonic_Rate)
            #print (WTS_Unique_Reads)
            #print ("HERE")
            #print (WTS_rRNA_contamination)
            #print (Yellow_Flags)
            #print (Red_Flags)

def parse_run_metrics(sample, run):
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics', run + '.align_bwa_mem.csv')
        with open(filename, 'r', encoding="utf-8") as file:
            for line in file:
                parsed_line = line.split(",")
                if sample == parsed_line[6]:
                    raw_reads_count = parsed_line[12]
                    raw_mean_coverage = parsed_line[41]
                    raw_median_insert_size = parsed_line[37]
                    raw_mean_insert_size = parsed_line[38]
                    raw_duplication_rate = parsed_line[15]
    except FileNotFoundError:
        raw_reads_count = None
        raw_mean_coverage = None
        raw_median_insert_size = None
        raw_mean_insert_size = None
        raw_duplication_rate = None
    return raw_reads_count, raw_mean_coverage, raw_median_insert_size, raw_mean_insert_size, raw_duplication_rate

def extract_rna_ribosomal(sample):
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample, 'rnaseqc', sample, sample + '.rRNA_counts.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            line = lines[0]
            fields = line.split("\t")
            ret = fields[0]
    except FileNotFoundError:
        ret = None
    return ret

# def extract_WTS_rRNA(ID):
#     if 'DT' in ID or 'DN' in ID or 'D-T' in ID or 'D-N' in ID:
#         return ''
#     else:
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna/' + ID + '/rnaseqc/*/*rRNA_counts.txt'
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 lines=f.readlines()
#                 line=lines[0]
#                 fields = line.split("\t")
#                 return fields[0]

def parse_rnaseqc_metrics_tmp(sample):
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample, 'rnaseqc', sample, sample + '.metrics.tmp.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            rna_aligned_reads_count = lines[3].split("\t")[0]
            rna_exonic_rate = round(float(lines[7].split("\t")[1])*100, 2)
    except FileNotFoundError:
        rna_aligned_reads_count = None
        rna_exonic_rate = None
    return rna_aligned_reads_count, rna_exonic_rate

# def extract_WTS_unique(ID):
#     if 'DT' in ID or 'DN' in ID or 'D-T' in ID or 'D-N' in ID:
#         return ''
#     else:
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna/' + ID + '/rnaseqc/*/*.metrics.tmp.txt'
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 lines=f.readlines()
#                 line=lines[3]
#                 fields = line.split("\t")
#                 Output = fields[0]
#                 return Output

# def extract_WTS_exonic(ID):
#     if 'DT' in ID or 'DN' in ID or 'D-T' in ID or 'D-N' in ID:
#         return ''
#     else:
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna/' + ID + '/rnaseqc/*/*.metrics.tmp.txt'
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 lines=f.readlines()
#                 line=lines[7]
#                 fields = line.split("\t")
#                 Output = float(fields[1])*100
#                 return (f"%.2f" % round(Output, 2))

# def extract_WTS_Clusters(ID):
#     if 'DT' in ID or 'DN' in ID or 'D-T' in ID or 'D-N' in ID:
#         return ''
#     else:
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/*'
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 for line in f:
#                     if ID in line:
#                         data = line.split(",")
#                         return data[12]

# def extract_raw_coverage(ID):
#     if 'R' in ID:
#         return ''
#     else:
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/*'
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 for line in f:
#                     if ID in line:
#                         data = line.split(",")
#                         return data[41]
# def extract_med_map_ins_size(ID):
#     path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/*'
#     for filename in glob.glob(path):
#         with open(filename, 'r') as f:
#             for line in f:
#                 if ID in line:
#                     data = line.split(",")
#                     return data[38]
# def extract_mean_map_ins_size(ID):
#     path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/*'
#     for filename in glob.glob(path):
#         with open(filename, 'r') as f:
#             for line in f:
#                 if ID in line:
#                     data = line.split(",")
#                     return data[39]

def extract_purity(sample, patient):
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants', patient, 'purple', sample + '.purple.purity.tsv')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            line = lines[1]
            fields = line.split("\t")
            ret = float(fields[0])*100
    except FileNotFoundError:
        ret = None
    return ret

# def extract_purity(ID,PATIENT):
#     if 'R' in ID or ID.endswith('N'):
#         return ''
#     else:
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/' + PATIENT + '/purple/*.purity.tsv'
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 lines=f.readlines()
#                 line=lines[1]
#                 fields = line.split("\t")
#                 Output = float(fields[0])*100
#                 return Output


def extract_contamination(patient, sample_type):
    if sample_type in ('DN', 'DT'):
        filename = "".join(os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics', patient + '.contamination.tsv'))
        try:
            with open(filename, 'r', encoding="utf-8") as file:
                for line in file:
                    if line.startswith('Normal') and sample_type == 'DN':
                        ret = line.split(" ")[-1][:-2]
                    elif line.startswith('Tumor') and sample_type == 'DT':
                        ret = line.split(" ")[-1][:-2]
        except FileNotFoundError:
            ret = None
    else:
        ret = None
    return ret

#WGS N % Contamination,WGS T % Contamination,
# def extract_contamination(ID,PATIENT):
#     if 'R' in ID:
#         return ''
#     else:
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/' + PATIENT + '*.contamination.tsv'
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 for line in f:
#                     words = line.split(" ")
#                     output = words[-1]
#                     output = output[:-2]
#                     if words[0].startswith('N') and ID.endswith('N'):
#                         return output
#                     elif words[0].startswith('T') and ID.endswith('T'):
#                         return output

def extract_concordance(patient, sample_type):
    if sample_type in ('DN', 'DT'):
        filename = "".join(os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics', patient + '.concordance.tsv'))
        try:
            with open(filename, 'r', encoding="utf-8") as file:
                for line in file:
                    if line.startswith('Concordance'):
                        ret = line.split(" ")[-1][:-2]
        except FileNotFoundError:
            ret = None
    else:
        ret = None
    return ret

# def extract_concordance(ID,PATIENT):
#     if 'R' in ID:
#         return ''
#     else:
#         CON = None
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/' + PATIENT + '*.concordance.tsv'
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 CON = f.readline()
#         if CON:
#             CON = CON.replace('Concordance:','')
#             CON = CON.replace(' ','')
#             CON = CON[:-1]
#             CON = CON[:-1]
#             if float(CON)<=1:
#                 CON=float(CON)*100
#             return CON
#         else:
#             return None


# def extract_sambama_dups(ID,PATIENT):
#     if 'R' in ID:
#         return ''
#     else:
#         DUPS = None;
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/run_metrics/*'
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 for line in f:
#                     if ID in line:
#                         data = line.split(",")
#                         DUPS = data[15]
#         if DUPS == None or DUPS == '':
#             return ''
#         else:
#             return float(DUPS)

def extract_insert_size(sample, patient, sample_type):
    try:
        if sample_type == 'RT':
            filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna', sample + '.insert_size_metrics')
            with open(filename, 'r', encoding="utf-8") as file:
                lines = file.readlines()
                line = lines[7]
                metrics = line.split("\t")
                ret = metrics[0]
        elif sample_type in ('DN', 'DT'):
            filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', patient + '.multiqc_data', 'multiqc_qualimap_bamqc_genome_results.txt')
            with open(filename, 'r', encoding="utf-8") as file:
                for line in file:
                    parsed_line = line.split("\t")
                    if parsed_line[0] == sample:
                        ret = round(float(parsed_line[7]), 0)
        else:
            ret = None
    except FileNotFoundError:
        ret = None
    return ret

# def extract_insert_size(ID,PATIENT):
#     if 'R' in ID:
#         OUTPUT = 0;
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna/' +  ID + '/*.insert_size_metrics'
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 lines=f.readlines()
#                 line=lines[7]
#                 metrics = line.split("\t")
#                 OUTPUT = metrics[0]
#         return OUTPUT
#     else:
#         return parse_multiqc2(ID,'median_insert_size',2,PATIENT)

def extract_dedup_coverage(sample):
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', sample, 'qualimap', sample, 'genome_results.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            lines = file.readlines()
            line = lines[71]
            #line is      mean coverageData = 304.9902X
            metrics = line.split(" ")
            ret = float(metrics[-1].replace('X', ''))
    except (FileNotFoundError, ValueError):
        ret = None
    return ret

# def extract_ded_coverage(ID,PATIENT):
#     if 'R' in ID:
#         return ''
#     else:
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna/' +  ID + '/qualimap/' + ID + '/genome_results.txt'
#         OUTPUT = ''
#         for filename in glob.glob(path):
#             with open(filename, 'r') as f:
#                 lines=f.readlines()
#                 line=lines[71]
#                 #line is      mean coverageData = 304.9902X
#                 metrics = line.split(" ")
#                 OUTPUT = metrics[-1].replace('X', '')
#         if OUTPUT == '' or OUTPUT == '':
#             return ""
#         else:
#             return float(OUTPUT)

def extract_min_aln_rds(sample, patient):
    ret = None
    try:
        filename = os.path.join('/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna', patient + '.multiqc_data', 'multiqc_general_stats.txt')
        with open(filename, 'r', encoding="utf-8") as file:
            for line in file:
                parsed_line = line.split("\t")
                if parsed_line[0] == sample:
                    ret = round(float(parsed_line[3]), 0)
    except FileNotFoundError:
        ret = None
    return ret

# def extract_min_aln_rds(ID,PATIENT):
#     if 'R' in ID:
#         return ''
#     else:
#         return parse_multiqc(ID,'QualiMap_mqc-generalstats-qualimap-mapped_reads',0,PATIENT)

# def parse_multiqc(ID,field,Round,PATIENT):
#     path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna/' +  PATIENT + '*.multiqc_data/multiqc_general_stats.txt'
#     for filename in glob.glob(path):
#         with open(filename, 'r') as f:
#             Output =''
#             header = f.readline().split("\t")
#             index = header.index(field)
#             A = f.readline().split("\t")
#             B = f.readline().split("\t")
#             if A[0].endswith(ID[-1]):
#                 Output = A[index]
#             elif B[0].endswith(ID[-1]):
#                 Output = B[index]
#             if Output != '':
#                 Output=float(Output)
#                 Output=f"%.{Round}f" % round(Output, Round)
#             return Output

# def parse_multiqc2(ID,field,Round,PATIENT):
#     path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna/' +  PATIENT + '*.multiqc_data/multiqc_qualimap_bamqc_genome_results.txt'
#     for filename in glob.glob(path):
#         with open(filename, 'r') as f:
#             Output =''
#             header = f.readline().split("\t")
#             index = header.index(field)
#             A = f.readline().split("\t")
#             B = f.readline().split("\t")
#             if A[0].endswith(ID[-1]):
#                 Output = A[index]
#             elif B[0].endswith(ID[-1]):
#                 Output = B[index]
#             if Output != '':
#                 Output=float(Output)
#                 Output=f"%.{Round}f" % round(Output, Round)
#             return Output

def extract_bs_over_q30(sample, sample_type):
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
        percent_abv = None
    return percent_abv

# def extract_bs_over_q30(ID,PATIENT):
#     path = ''
#     if 'R' in ID:
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna/' + ID + '/picard_metrics/*quality_distribution_metrics'
#     else:
#         path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna/' + ID + '*/picard_metrics/*quality_distribution_metrics'
#     Tester = re.compile('(\d+)\W+(\d+)')
#     for filename in glob.glob(path):
#         with open(filename, 'r') as f:
#             Abv_30=0
#             Blw_30=0
#             for line in f:
#                 if line[:1].isdigit():
#                     Test = Tester.match(line)
#                     Qual = Test.group(1)
#                     count = Test.group(2)
#                     if (int(Qual) < 30):
#                         Blw_30 += int(count)
#                     else :
#                         Abv_30 += int(count)
#             percent_abv=(Abv_30/(Abv_30+Blw_30))*100
#             percent_abv = "%.2f" % round(percent_abv, 2)
#             return percent_abv


if __name__ == '__main__':
    main()
    #Update db with the objects
