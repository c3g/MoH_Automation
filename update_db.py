#!/usr/bin/env python3

import progressbar
import  moh_resources

WIDGETS = [' [', progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') - ', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']

def main():
    beluga_moh_folder = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING"
    # beluga_db = os.path.join(beluga_moh_folder, "DATABASE", MOH_analysis.db)
    # TEST
    beluga_db = "/Users/pstretenowich/Documents/local/projects/MOH/MOH_analysis.db"
    connection = moh_resources.create_connection(beluga_db)

    patients = moh_resources.extract_patient_names(connection)
    samples_list = []
    paired_samples_dict = {}

    print("Updating Database by patient...")
    with progressbar.ProgressBar(max_value=len(patients), widgets=WIDGETS) as progress:
        for index, patient in enumerate(patients, 1):
            run_sample = moh_resources.extract_patient_field(connection, patient, "run")
            dna_n_sample = moh_resources.extract_patient_field(connection, patient, "dna_n")
            if dna_n_sample:
                samples_list.append(dna_n_sample)
                paired_samples_dict[dna_n_sample] = (patient, run_sample)
            dna_t_sample = moh_resources.extract_patient_field(connection, patient, "dna_t")
            if dna_t_sample:
                samples_list.append(dna_t_sample)
                paired_samples_dict[dna_t_sample] = (patient, run_sample)
            rna_sample = moh_resources.extract_patient_field(connection, patient, "rna")
            if rna_sample:
                samples_list.append(rna_sample)
                paired_samples_dict[rna_sample] = (patient, run_sample)
            # Creating patient object
            patient_data = moh_resources.Progress(connection, patient, beluga_moh_folder)
            # Updating timestamp table
            moh_resources.update_timestamp_details(patient_data)
            # Updating file_location table
            moh_resources.update_fileloc_details(patient_data)
            patient_data.gather_bam_loc()
            patient_data.gather_final_bams()
            patient_data.gather_vcfs()
            patient_data.gather_reports()
            patient_data.gather_rna_other()
            moh_resources.update_timestamp_details(patient_data)
            moh_resources.update_fileloc_details(patient_data)
            # Updating status table
            patient_data.update_status()
            progress.update(index)
    extract_data(samples_list, connection, paired_samples_dict)
    print("Committing changes to Database...")
    connection.commit()
    connection.close()
    print("Done.")

def extract_data(samples_list, connection, paired_samples_dict):
    print("Updating Database by sample...")
    with progressbar.ProgressBar(max_value=len(samples_list), widgets=WIDGETS) as progress:
        for index, sample in enumerate(samples_list, 1):
            print(sample + "\n")
            patient = paired_samples_dict[sample][0]
            run = paired_samples_dict[sample][1]
            if sample.endswith('DN'):
                sample_type = 'DN'
            elif sample.endswith('DT'):
                sample_type = 'DT'
            elif sample.endswith('RT'):
                sample_type = 'RT'
            flag = []
            fail = []

            dna_bases_over_q30_percent = moh_resources.extract_bs_over_q30(sample, sample_type)
            moh_resources.dna_bases_over_q30_percent_check(dna_bases_over_q30_percent, sample_type, fail, flag)

            dna_aligned_reads_count = moh_resources.extract_min_aln_rds(sample, patient)
            moh_resources.dna_aligned_reads_count_check(dna_aligned_reads_count, sample_type, fail, flag)

            dna_dedup_coverage = moh_resources.extract_dedup_coverage(sample)

            median_insert_size = moh_resources.extract_insert_size(sample, patient, sample_type)
            moh_resources.median_insert_size_check(median_insert_size, fail, flag)

            dna_contamination = moh_resources.extract_contamination(patient, sample_type)
            moh_resources.dna_contamination_check(dna_contamination, fail)

            dna_concordance = moh_resources.extract_concordance(patient, sample_type)
            moh_resources.dna_concordance_check(dna_concordance, fail)

            dna_tumour_purity = moh_resources.extract_purity(sample, patient)
            moh_resources.dna_tumour_purity_check(dna_tumour_purity, fail)

            rna_aligned_reads_count, rna_exonic_rate = moh_resources.parse_rnaseqc_metrics_tmp(sample)
            moh_resources.rna_exonic_rate_check(rna_exonic_rate, fail, flag)

            rrna_count = moh_resources.extract_rna_ribosomal(sample)
            rna_ribosomal_contamination_count = moh_resources.rna_ribosomal_contamination_count_check(rrna_count, rna_aligned_reads_count, fail, flag)

            # Flags
            flags.extend(moh_resources.extract_value(connection, "key_metric", sample, "flag").split(";"))
            fails.extend(moh_resources.extract_value(connection, "key_metric", sample, "fail").split(";"))
            flags = ';'.join(set(flag))
            fails = ';'.join(set(fail))

            moh_resources.update_key_metric_table(
                conn=connection,
                sample=sample,
                dna_bases_over_q30_percent=dna_bases_over_q30_percent,
                dna_aligned_reads_count=dna_aligned_reads_count,
                dna_dedup_coverage=dna_dedup_coverage,
                median_insert_size=median_insert_size,
                dna_contamination=dna_contamination,
                rna_aligned_reads_count=rna_aligned_reads_count,
                rna_exonic_rate=rna_exonic_rate,
                rna_ribosomal_contamination_count=rna_ribosomal_contamination_count,
                dna_concordance=dna_concordance,
                dna_tumour_purity=dna_tumour_purity,
                flag=flags,
                fail=fails
                )
            progress.update(index)

if __name__ == '__main__':
    main()
    #Update db with the objects
