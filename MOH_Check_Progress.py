#!/usr/bin/env python3
import progressbar
import glob
import sys
import re
import os.path
import time
import datetime
import errno
import logging
from  DB_OPS import update_metrics_db,create_connection,extract_sample_details,extract_fileloc_details,extract_timestamp_details,update_timestamp_details,update_fileloc_details,extract_sample_names,update_status_db


WIDGETS = [' [', progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') - ', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']

logger = logging.getLogger(__name__)

class SampleData:
    def __init__(self, connection, sample):
        self.conn = connection
        sample_details = extract_sample_details(connection, sample)
        self.sample = sample_details["Sample"]
        self.sample_true = sample_details["Sample_True"]
        self.institution = sample_details["Instituion"]
        self.cohort = sample_details["Cohort"]
        self.dna_n = sample_details["DNA_N"]
        self.dna_n_true = sample_details["DNA_N_True"]
        self.dna_t = sample_details["DNA_T"]
        self.dna_t_true = sample_details["DNA_T_True"]
        self.rna = sample_details["RNA"]
        self.rna_true = sample_details["RNA_True"]
    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

class Progress(SampleData):
    def grab_db_values(self):
        data = extract_fileloc_details(self.conn, self.sample)
        self.run_proc_bam_dna_t = data["Run_Proc_BAM_DNA_T"]
        self.run_proc_bam_dna_n = data["Run_Proc_BAM_DNA_N"]
        self.beluga_bam_dna_t = data["Beluga_BAM_DNA_T"]
        self.beluga_bam_dna_n = data["Beluga_BAM_DNA_N"]
        self.dna_vcf_g = data["DNA_VCF_G"]
        self.dna_vcf_s = data["DNA_VCF_S"]
        self.mutect2_somatic_vcf = data["Mutect2_Somatic_vcf"]
        self.strelka2_germline_vcf = data["strelka2_Germline_vcf"]
        self.strelka2_somatic_vcf = data["strelka2_Somatic_vcf"]
        self.vardict_germline_vcf = data["vardict_Germline_vcf"]
        self.vardict_somatic_vcf = data["vardict_Somatic_vcf"]
        self.varscan2_germline_vcf = data["varscan2_Germline_vcf"]
        self.varscan2_somatic_vcf = data["varscan2_Somatic_vcf"]
        self.gripss_somatic = data["gripss_somatic"]
        self.gripss_germline = data["gripss_germline"]
        self.purple_somatic = data["purple_somatic"]
        self.purple_germline = data["purple_germline"]
        self.purple_circos = data["purple_circos"]
        self.cnvkit_vcf = data["cnvkit_vcf"]
        self.final_vcf = data["Final_VCF"]
        self.final_dna_bam_t = data["Final_DNA_BAM_T"]
        self.final_dna_bam_n = data["Final_DNA_BAM_N"]
        self.dna_multiqc = data["DNA_MultiQC"]
        self.pcgr_report = data["pcgr_report"]
        self.pcgr_maf = data["pcgr_maf"]
        self.pcgr_snvs_indels = data["pcgr_snvs_indels"]
        self.pcgr_cna_segments = data["pcgr_cna_segments"]
        self.tp_ini = data["TP_ini"]
        self.run_proc_fastq_1_rna = data["Run_Proc_fastq_1_RNA"]
        self.run_proc_fastq_2_rna = data["Run_Proc_fastq_2_RNA"]
        self.beluga_fastq_1_rna = data["Beluga_fastq_1_RNA"]
        self.beluga_fastq_2_rna = data["Beluga_fastq_2_RNA"]
        self.rna_vcf = data["RNA_VCF"]
        self.rna_vcf_filt = data["RNA_VCF_filt"]
        self.final_rna_bam = data["Final_RNA_BAM"]
        self.rna_multiqc = data["RNA_MultiQC"]
        self.rna_pcgr_report = data["rna_pcgr_report"]
        self.rna_pcgr_maf = data["rna_pcgr_maf"]
        self.rna_pcgr_snvs_indels = data["rna_pcgr_snvs_indels"]
        self.rna_pcgr_cna_segments = data["rna_pcgr_cna_segments"]
        self.annofuse = data["AnnoFuse"]
        self.gridss = data["GRIDSS"]
        self.rna_abundance = data["RNA_Abundance"]
        self.big_wig_tracks_f = data["big_wig_tracks_F"]
        self.big_wig_tracks_r = data["big_wig_tracks_R"]
        self.rna_abundance_ini = data["RNA_Abundance_ini"]
        self.rna_variants_ini = data["RNA_Variants_ini"]
        #get the old timestamps
        data = extract_timestamp_details(self.conn, self.sample)
        self.ts_run_proc_bam_dna_t = data["Run_Proc_BAM_DNA_T"]
        self.ts_run_proc_bam_dna_n = data["Run_Proc_BAM_DNA_N"]
        self.ts_beluga_bam_dna_t = data["Beluga_BAM_DNA_T"]
        self.ts_beluga_bam_dna_n = data["Beluga_BAM_DNA_N"]
        self.ts_dna_vcf_g = data["DNA_VCF_G"]
        self.ts_dna_vcf_s = data["DNA_VCF_S"]
        self.ts_mutect2_somatic_vcf = data["Mutect2_Somatic_vcf"]
        self.ts_strelka2_germline_vcf = data["strelka2_Germline_vcf"]
        self.ts_strelka2_somatic_vcf = data["strelka2_Somatic_vcf"]
        self.ts_vardict_germline_vcf = data["vardict_Germline_vcf"]
        self.ts_vardict_somatic_vcf = data["vardict_Somatic_vcf"]
        self.ts_varscan2_germline_vcf = data["varscan2_Germline_vcf"]
        self.ts_varscan2_somatic_vcf = data["varscan2_Somatic_vcf"]
        self.ts_gripss_somatic = data["gripss_somatic"]
        self.ts_gripss_germline = data["gripss_germline"]
        self.ts_purple_somatic = data["purple_somatic"]
        self.ts_purple_germline = data["purple_germline"]
        self.ts_purple_circos = data["purple_circos"]
        self.ts_cnvkit_vcf = data["cnvkit_vcf"]
        self.ts_final_vcf = data["Final_VCF"]
        self.ts_final_dna_bam_t = data["Final_DNA_BAM_T"]
        self.ts_final_dna_bam_n = data["Final_DNA_BAM_N"]
        self.ts_dna_multiqc = data["DNA_MultiQC"]
        self.ts_pcgr_report = data["pcgr_report"]
        self.ts_pcgr_maf = data["pcgr_maf"]
        self.ts_pcgr_snvs_indels = data["pcgr_snvs_indels"]
        self.ts_pcgr_cna_segments = data["pcgr_cna_segments"]
        self.ts_tp_ini = data["TP_ini"]
        self.ts_run_proc_fastq_1_rna = data["Run_Proc_fastq_1_RNA"]
        self.ts_run_proc_fastq_2_rna = data["Run_Proc_fastq_2_RNA"]
        self.ts_beluga_fastq_1_rna = data["Beluga_fastq_1_RNA"]
        self.ts_beluga_fastq_2_rna = data["Beluga_fastq_2_RNA"]
        self.ts_rna_vcf = data["RNA_VCF"]
        self.ts_rna_vcf_filt = data["RNA_VCF_filt"]
        self.ts_final_rna_bam = data["Final_RNA_BAM"]
        self.ts_rna_multiqc = data["RNA_MultiQC"]
        self.ts_rna_pcgr_report = data["rna_pcgr_report"]
        self.ts_rna_pcgr_maf = data["rna_pcgr_maf"]
        self.ts_rna_pcgr_snvs_indels = data["rna_pcgr_snvs_indels"]
        self.ts_rna_pcgr_cna_segments = data["rna_pcgr_cna_segments"]
        self.ts_annofuse = data["AnnoFuse"]
        self.ts_gridss = data["GRIDSS"]
        self.ts_rna_abundance = data["RNA_Abundance"]
        self.ts_big_wig_tracks_f = data["big_wig_tracks_F"]
        self.ts_big_wig_tracks_r = data["big_wig_tracks_R"]
        self.ts_rna_abundance_ini = data["RNA_Abundance_ini"]
        self.ts_rna_variants_ini = data["RNA_Variants_ini"]

    def __init__(self, connection, sample):
        super().__init__(connection, sample)
        self.grab_db_values()

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

#################ONCE DONE JUST STORE THE WHOLE FLIPPEN INI HERE ########################################
    def Gather_DNA_ini(self):
        filename = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/TumorPair.config.trace.ini'
        self.tp_ini = filename
        self.ts_tp_ini = getime(filename)
#########################################################################################################
#################ONCE DONE JUST STORE THE WHOLE FLIPPEN INI HERE########################################
    def Gather_RNA_Light_ini(self):
        rna_abundance_ini = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/RnaSeqLight.config.trace.ini'
        self.rna_abundance_ini = rna_abundance_ini
        self.ts_rna_abundance_ini = getime(rna_abundance_ini)

    def Gather_RNA_Variants_ini(self):
        rna_variants_ini = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/RnaSeq.config.trace.ini'
        self.rna_variants_ini = rna_variants_ini
        self.ts_rna_variants_ini = getime(rna_variants_ini)
#########################################################################################################

    def Update_status(self):
        # DNA
        dna_n_transferred = "NA"
        if self.ts_run_proc_bam_dna_n != "NA":
            dna_n_transferred = "Complete"
        dna_t_transferred = "NA"
        if self.ts_run_proc_bam_dna_t != "NA":
            dna_t_transferred = "Complete"
        dna_alignment= "NA"
        if self.final_dna_bam_t != "NA":
            dna_alignment = "Complete"
        dna_variant_call = "NA"
        if self.dna_vcf_g != "NA" and self.dna_vcf_s != "NA" :
            dna_variant_call = "Complete"
        dna_report = "NA"
        if self.dna_multiqc != "NA" and self.pcgr_report != "NA":
            dna_report = "Complete"
        dna_pipeline_execution = "NA"
        if dna_n_transferred == dna_t_transferred == dna_alignment == dna_variant_call == dna_report == "Complete":
            dna_pipeline_execution = "Complete"
        # RNA
        rna_transferred = "NA"
        if self.beluga_fastq_1_rna != "NA":
            rna_transferred = "Complete"
        # For abundance
        rna_pipeline_light_execution = "NA"
        if self.rna_abundance != "NA":
            rna_pipeline_light_execution = "Complete"
        rna_alignment = "NA"
        if self.final_rna_bam != "NA":
            rna_alignment = "Complete"
        rna_variant_call = "NA"
        # if self.rna_vcf != "NA" and self.rna_vcf_filt != "NA":
        if self.rna_vcf != "NA":
            rna_variant_call = "Complete"
        rna_report = "NA"
        if self.annofuse != "NA" and self.rna_multiqc != "NA" and self.rna_pcgr_report != "NA":
            rna_report = "Complete"
        rna_pipeline_execution = "NA"
        if rna_transferred == rna_alignment == rna_variant_call == rna_report == "Complete":
            rna_pipeline_execution = "Complete"
        # DELIVERY
        dna_delivered, rna_light_delivered, rna_delivered = self.check_delivery()
        update_status_db(
            self.conn,
            self.sample,
            dna_n_transferred,
            dna_t_transferred,
            dna_alignment,
            dna_variant_call,
            dna_report,
            dna_pipeline_execution,
            dna_delivered,
            rna_transferred,
            rna_pipeline_light_execution,
            rna_light_delivered,
            rna_alignment,
            rna_variant_call,
            rna_report,
            rna_pipeline_execution,
            rna_delivered
            )

    def check_delivery(self):
        """
        Checks delivery status for Tumour Pair Pipeline, RNASeq Light Pipeline and RNASeq Variants Pipeline
        """
        dna_delivered = rna_light_delivered = rna_delivered = "NA"
        base_dir = "/lustre03/project/6007512/C3G/projects/share/MOH"
        patient_dir = os.path.join(base_dir, self.institution, self.cohort, self.sample, "parameters")
        if os.path.isfile(os.path.join(patient_dir, f"{self.sample}.TumourPair.ini")):
            dna_delivered = "Complete"
        if os.path.isfile(os.path.join(patient_dir, f"{self.sample}.RNA.Light.ini")):
            rna_light_delivered = "Complete"
        if os.path.isfile(os.path.join(patient_dir, f"{self.sample}.RNA.Variants.ini")):
            rna_delivered = "Complete"
        return dna_delivered, rna_light_delivered, rna_delivered

    def Gather_Run_Proc_BAM(self):
        if self.dna_n_true == "NA":
            self.run_proc_bam_dna_t = "NA"
            self.run_proc_bam_dna_n = "NA"
            self.ts_run_proc_bam_dna_t = "NA"
            self.ts_run_proc_bam_dna_n = "NA"
        else:
            path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*'
            dna_n = False
            dna_t = False
            if self.run_proc_bam_dna_t == "NA" or self.run_proc_bam_dna_n == "NA":
                for filename in glob.glob(path):
                    with open(filename, 'r') as file:
                        for line in file:
                            if dna_n and dna_t:
                                break
                            if self.dna_n in line:
                                fields = line.split(",")
                                self.run_proc_bam_dna_n = fields[0].strip()
                                self.ts_run_proc_bam_dna_n  = getime(filename)
                                dna_n = True
                            elif self.dna_t in line:
                                fields = line.split(",")
                                self.run_proc_bam_dna_t = fields[0].strip()
                                self.ts_run_proc_bam_dna_t = getime(filename)
                                dna_t = True
                if not dna_n:
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"Noting found for {dna_n} in /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*")
                if not dna_t:
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"Noting found for {dna_t} in /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*")

        if self.rna_true == "NA":
            self.run_proc_fastq_1_rna = "NA"
            self.run_proc_fastq_2_rna = "NA"
            self.ts_run_proc_fastq_1_rna = "NA"
            self.ts_run_proc_fastq_2_rna = "NA"
        else:
            path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*'
            fastq1 = False
            fastq2 = False
            if self.run_proc_fastq_1_rna == "NA" or self.run_proc_fastq_2_rna == "NA":
                for filename in glob.glob(path):
                    with open(filename, 'r') as file:
                        for line in file:
                            if fastq1 and fastq2:
                                break
                            if self.rna in line:
                                fields = line.split(",")
                                if re.search(r"R1.*.fastq", os.path.basename(fields[0])):
                                    self.run_proc_fastq_1_rna = fields[0].strip()
                                    self.ts_run_proc_fastq_1_rna = getime(filename)
                                    fastq1 = True
                                elif re.search(r"R2.*.fastq", os.path.basename(fields[0])):
                                    self.run_proc_fastq_2_rna = fields[0].strip()
                                    self.ts_run_proc_fastq_2_rna = getime(filename)
                                    fastq2 = True
                if not fastq1:
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"Noting found for {fastq1} in /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*")
                if not fastq2:
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"Noting found for {fastq2} in /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*")

    def Gather_BAM_loc(self):
        loc1 = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/raw_reads"
        loc2 = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/raw_reads"
        if self.dna_n_true == "NA" or self.dna_t_true == "NA":
            self.beluga_bam_dna_t = "NA"
            self.beluga_bam_dna_n= "NA"
            self.ts_beluga_bam_dna_t = "NA"
            self.ts_beluga_bam_dna_n = "NA"
        else:
            # bam dna_n in raw_reads
            dna_n_file_rr = os.path.join(loc1, self.dna_n, os.path.basename(self.run_proc_bam_dna_n))
            dna_n_file_old_rr = os.path.join(loc1, self.dna_n, self.dna_n + ".bam")
            if os.path.exists(dna_n_file_rr):
                self.beluga_bam_dna_n = dna_n_file_rr
                self.ts_beluga_bam_dna_n = getime(dna_n_file_rr)
            elif os.path.exists(dna_n_file_old_rr):
                self.beluga_bam_dna_n = dna_n_file_old_rr
                self.ts_beluga_bam_dna_n = getime(dna_n_file_old_rr)
            # bam dna_t in raw_reads
            dna_t_file_rr = os.path.join(loc1, self.dna_t, os.path.basename(self.run_proc_bam_dna_t))
            dna_t_file_old_rr = os.path.join(loc1, self.dna_t, self.dna_t + ".bam")
            if os.path.exists(dna_t_file_rr):
                self.beluga_bam_dna_t = dna_t_file_rr
                self.ts_beluga_bam_dna_t = getime(dna_t_file_rr)
            elif os.path.exists(dna_t_file_old_rr):
                self.beluga_bam_dna_t = dna_t_file_old_rr
                self.ts_beluga_bam_dna_t = getime(dna_t_file_old_rr)
            # bam dna_n in MAIN/raw_reads
            dna_n_file = os.path.join(loc2, self.dna_n, os.path.basename(self.run_proc_bam_dna_n))
            dna_n_file_old = os.path.join(loc2, self.dna_n, self.dna_n + ".bam")
            if os.path.exists(dna_n_file):
                self.beluga_bam_dna_n = dna_n_file
                self.ts_beluga_bam_dna_n = getime(dna_n_file)
            elif os.path.exists(dna_n_file_old):
                self.beluga_bam_dna_n = dna_n_file_old
                self.ts_beluga_bam_dna_n = getime(dna_n_file_old)
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"{dna_n_file_rr} nor {dna_n_file_old_rr} nor {dna_n_file} nor {dna_n_file_old}")
            # bam dna_t in MAIN/raw_reads
            dna_t_file = os.path.join(loc2, self.dna_t, os.path.basename(self.run_proc_bam_dna_t))
            dna_t_file_old = os.path.join(loc2, self.dna_t, self.dna_t + ".bam")
            if os.path.exists(dna_t_file):
                self.beluga_bam_dna_t = dna_t_file
                self.ts_beluga_bam_dna_t = getime(dna_t_file)
            elif os.path.exists(dna_t_file_old):
                self.beluga_bam_dna_t = dna_t_file_old
                self.ts_beluga_bam_dna_t = getime(dna_t_file_old)
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"{dna_t_file_rr} nor {dna_t_file_old_rr} nor {dna_t_file} nor {dna_t_file_old}")
        if self.rna_true == "NA":
            self.beluga_fastq_1_rna = "NA"
            self.beluga_fastq_2_rna = "NA"
            self.ts_beluga_fastq_1_rna = "NA"
            self.ts_beluga_fastq_2_rna = "NA"
        else:
            # fastq1 in raw_reads
            rna_fq1_rr = os.path.join(loc1, self.rna, os.path.basename(self.run_proc_fastq_1_rna))
            rna_fq1_old_rr = os.path.join(loc1, self.rna, self.rna + "_R1.fastq.gz")
            if os.path.exists(rna_fq1_rr):
                self.beluga_fastq_1_rna = rna_fq1_rr
                self.ts_beluga_fastq_1_rna = getime(rna_fq1_rr)
            elif os.path.exists(rna_fq1_old_rr):
                self.beluga_fastq_1_rna = rna_fq1_old_rr
                self.ts_beluga_fastq_1_rna = getime(rna_fq1_old_rr)
            # fastq2 in raw_reads
            rna_fq2_rr = os.path.join(loc1, self.rna, os.path.basename(self.run_proc_fastq_2_rna))
            rna_fq2_old_rr = os.path.join(loc1, self.rna, self.rna + "_R2.fastq.gz")
            if os.path.exists(rna_fq2_rr):
                self.beluga_fastq_2_rna = rna_fq2_rr
                self.ts_beluga_fastq_2_rna = getime(rna_fq2_rr)
            elif os.path.exists(rna_fq2_old_rr):
                self.beluga_fastq_2_rna = rna_fq2_old_rr
                self.ts_beluga_fastq_2_rna = getime(rna_fq2_old_rr)
            # fastq1 in MAIN/raw_reads
            rna_fq1 = os.path.join(loc2, self.rna, os.path.basename(self.run_proc_fastq_1_rna))
            rna_fq1_old = os.path.join(loc2, self.rna, self.rna + "_R1.fastq.gz")
            if os.path.exists(rna_fq1):
                self.beluga_fastq_1_rna = rna_fq1
                self.ts_beluga_fastq_1_rna = getime(rna_fq1)
            elif os.path.exists(rna_fq1_old):
                self.beluga_fastq_1_rna = rna_fq1_old
                self.ts_beluga_fastq_1_rna = getime(rna_fq1_old)
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"{rna_fq1_rr} nor {rna_fq1_old_rr} nor {rna_fq1_old} nor {rna_fq1_old}")
            # fastq2 in MAIN/raw_reads
            rna_fq2 = os.path.join(loc2, self.rna, os.path.basename(self.run_proc_fastq_2_rna))
            rna_fq2_old = os.path.join(loc2, self.rna, self.rna + "_R2.fastq.gz")
            if os.path.exists(rna_fq2):
                self.beluga_fastq_2_rna = rna_fq2
                self.ts_beluga_fastq_2_rna = getime(rna_fq2)
            elif os.path.exists(rna_fq2_old):
                self.beluga_fastq_2_rna = rna_fq2_old
                self.ts_beluga_fastq_2_rna = getime(rna_fq2_old)
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"{rna_fq2_rr} nor {rna_fq2_old_rr} nor {rna_fq2_old} nor {rna_fq2_old}")

    def Gather_Final_BAMs(self):
        if self.dna_n_true == "NA":
            self.final_dna_bam_t = "NA"
            self.final_dna_bam_n = "NA"
            self.ts_final_dna_bam_t= "NA"
            self.ts_final_dna_bam_n= "NA"
        else:
            final_dna_bam_n_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.dna_n, self.dna_n + ".sorted.dup.recal.bam")
            try:
                if self.ts_final_dna_bam_n != getime(final_dna_bam_n_file):
                    self.final_dna_bam_n = final_dna_bam_n_file
                    self.ts_final_dna_bam_n = getime(final_dna_bam_n_file)
                    self.Gather_DNA_ini()
            except FileNotFoundError:
                pass
            final_dna_bam_t_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.dna_t, self.dna_t + ".sorted.dup.recal.bam")
            try:
                if self.ts_final_dna_bam_t != getime(final_dna_bam_t_file):
                    self.final_dna_bam_t = final_dna_bam_t_file
                    self.ts_final_dna_bam_t = getime(final_dna_bam_t_file)
                    self.Gather_DNA_ini()
            except FileNotFoundError:
                pass
        if self.rna_true == "NA":
            self.final_rna_bam = "NA"
            self.ts_final_rna_bam = "NA"
        else:
            final_rna_bam_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, self.rna + ".sorted.mdup.split.recal.bam")
            try:
                if self.ts_final_rna_bam != getime(final_rna_bam_file):
                    self.ts_final_rna_bam = getime(final_rna_bam_file)
                    self.final_rna_bam_variants = final_rna_bam_file
            except FileNotFoundError:
                pass

    def Gather_VCFs(self):
        if self.dna_n_true == "NA":
            self.dna_vcf_g = "NA"
            self.dna_vcf_s = "NA"
            self.mutect2_somatic_vcf = "NA"
            self.strelka2_germline_vcf = "NA"
            self.strelka2_somatic_vcf = "NA"
            self.vardict_germline_vcf = "NA"
            self.vardict_somatic_vcf = "NA"
            self.varscan2_germline_vcf = "NA"
            self.varscan2_somatic_vcf = "NA"
            self.cnvkit_vcf = "NA"
            self.final_vcf = "NA"
            self.ts_dna_vcf_g = "NA"
            self.ts_dna_vcf_s = "NA"
            self.ts_mutect2_somatic_vcf = "NA"
            self.ts_strelka2_germline_vcf = "NA"
            self.ts_strelka2_somatic_vcf = "NA"
            self.ts_vardict_germline_vcf = "NA"
            self.ts_vardict_somatic_vcf = "NA"
            self.ts_varscan2_germline_vcf = "NA"
            self.ts_varscan2_somatic_vcf = "NA"
            self.ts_cnvkit_vcf = "NA"
            self.ts_final_vcf = "NA"
        else:
            #####################Not Implemented########################################
            self.ts_final_vcf = "NA"
            self.final_vcf = "NA"
            ############################################################################
            dna_vcf_g_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble", self.sample, self.sample + ".ensemble.germline.vt.annot.vcf.gz")
            try:
                if self.ts_dna_vcf_g != getime(dna_vcf_g_file):
                    self.ts_dna_vcf_g = getime(dna_vcf_g_file)
                    self.dna_vcf_g = dna_vcf_g_file
            except FileNotFoundError:
                pass
            dna_vcf_s_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble", self.sample, self.sample + ".ensemble.somatic.vt.annot.vcf.gz")
            try:
                if self.ts_dna_vcf_s != getime(dna_vcf_s_file):
                    self.ts_dna_vcf_s = getime(dna_vcf_s_file)
                    self.dna_vcf_s = dna_vcf_s_file
            except FileNotFoundError:
                pass
            mutect2_somatic_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants", self.sample, self.sample + ".mutect2.somatic.vt.vcf.gz")
            try:
                if self.ts_mutect2_somatic_vcf != getime(mutect2_somatic_vcf_file):
                    self.ts_mutect2_somatic_vcf = getime(mutect2_somatic_vcf_file)
                    self.mutect2_somatic_vcf = mutect2_somatic_vcf_file
            except FileNotFoundError:
                pass
            strelka2_somatic_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants", self.sample, self.sample + ".strelka2.somatic.purple.vcf.gz")
            try:
                if self.ts_strelka2_somatic_vcf != getime(strelka2_somatic_vcf_file):
                    self.ts_strelka2_somatic_vcf = getime(strelka2_somatic_vcf_file)
                    self.strelka2_somatic_vcf = strelka2_somatic_vcf_file
            except FileNotFoundError:
                pass
            strelka2_germline_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants", self.sample, self.sample + ".strelka2.germline.vt.vcf.gz")
            try:
                if self.ts_strelka2_germline_vcf != getime(strelka2_germline_vcf_file):
                    self.ts_strelka2_germline_vcf = getime(strelka2_germline_vcf_file)
                    self.strelka2_germline_vcf = strelka2_germline_vcf_file
            except FileNotFoundError:
                pass
            vardict_germline_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants", self.sample, self.sample + ".vardict.germline.vt.vcf.gz")
            try:
                if self.ts_vardict_germline_vcf != getime(vardict_germline_vcf_file):
                    self.ts_vardict_germline_vcf = getime(vardict_germline_vcf_file)
                    self.vardict_germline_vcf= vardict_germline_vcf_file
            except FileNotFoundError:
                pass
            vardict_somatic_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants", self.sample, self.sample + ".vardict.somatic.vt.vcf.gz")
            try:
                if self.ts_vardict_somatic_vcf != getime(vardict_somatic_vcf_file):
                    self.ts_vardict_somatic_vcf = getime(vardict_somatic_vcf_file)
                    self.vardict_somatic_vcf = vardict_somatic_vcf_file
            except FileNotFoundError:
                pass
            varscan2_germline_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants", self.sample, self.sample + ".varscan2.germline.vt.vcf.gz")
            try:
                if self.ts_varscan2_germline_vcf != getime(varscan2_germline_vcf_file):
                    self.ts_varscan2_germline_vcf = getime(varscan2_germline_vcf_file)
                    self.varscan2_germline_vcf = varscan2_germline_vcf_file
            except FileNotFoundError:
                pass
            varscan2_somatic_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants", self.sample, self.sample + ".varscan2.somatic.vt.vcf.gz")
            try:
                if self.ts_varscan2_somatic_vcf != getime(varscan2_somatic_vcf_file):
                    self.ts_varscan2_somatic_vcf = getime(varscan2_somatic_vcf_file)
                    self.varscan2_somatic_vcf = varscan2_somatic_vcf_file
            except FileNotFoundError:
                pass
            cnvkit_vcf = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/SVariants", self.sample, self.sample + ".cnvkit.vcf.gz")
            try:
                if self.ts_cnvkit_vcf != getime(cnvkit_vcf):
                    self.ts_cnvkit_vcf = getime(cnvkit_vcf)
                    self.cnvkit_vcf = cnvkit_vcf
            except FileNotFoundError:
                pass
        if self.rna_true == "NA":
            self.rna_vcf = "NA"
            self.ts_rna_vcf = "NA"
            self.rna_vcf_filt = "NA"
            self.ts_rna_vcf_filt = "NA"
        else:
            rna_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, self.rna + ".hc.vt.annot.vcf.gz")
            try:
                if self.ts_rna_vcf != getime(rna_vcf_file):
                    self.rna_vcf = rna_vcf_file
                    self.ts_rna_vcf = getime(rna_vcf_file)
                    self.Gather_RNA_Variants_ini()
            except FileNotFoundError:
                pass
            rna_vcf_filt_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, self.rna + ".hc.vt.annot.flt.vcf.gz")
            try:
                if self.ts_rna_vcf_filt != getime(rna_vcf_filt_file):
                    self.rna_vcf_filt = rna_vcf_filt_file
                    self.ts_rna_vcf_filt = getime(rna_vcf_filt_file)
                    self.Gather_RNA_Variants_ini()
            except FileNotFoundError:
                pass


    def Gather_reports(self):
        if self.dna_n_true == "NA":
            self.dna_multiqc = "NA"
            self.ts_dna_multiqc = "NA"
            self.ts_pcgr = "NA"
        else:
            dna_multiqc_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna", self.sample + ".multiqc.html")
            try:
                if self.ts_dna_multiqc != getime(dna_multiqc_file):
                    self.ts_dna_multiqc = getime(dna_multiqc_file)
                    self.dna_multiqc = dna_multiqc_file
            except FileNotFoundError:
                pass

        if self.rna_true == "NA":
            self.rna_multiqc = "NA"
            self.ts_rna_multiqc = "NA"
        else:
            rna_multiqc_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna", self.sample + ".multiqc.html")
            try:
                if self.ts_rna_multiqc != getime(rna_multiqc_file):
                    self.ts_rna_multiqc = getime(rna_multiqc_file)
                    self.rna_multiqc = rna_multiqc_file
            except FileNotFoundError:
                pass
            self.rna_multiqc = "NA"
            self.ts_rna_multiqc = "NA"

    def Gather_PCGR(self):
        if self.dna_n_true == "NA":
            self.pcgr_report = "NA"
            self.pcgr_maf = "NA"
            self.pcgr_snvs_indels = "NA"
            self.pcgr_cna_segments = "NA"
            self.ts_pcgr_report = "NA"
            self.ts_pcgr_maf = "NA"
            self.ts_pcgr_snvs_indels = "NA"
            self.ts_pcgr_cna_segments = "NA"
        else:
            pcgr_report = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble", self.sample, "pcgr", self.sample + ".pcgr_acmg.grch38.flexdb.html")
            try:
                if self.ts_pcgr_report != getime(pcgr_report):
                    self.ts_pcgr_report = getime(pcgr_report)
                    self.pcgr_report = pcgr_report
            except FileNotFoundError:
                pass
            pcgr_maf = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble", self.sample, "pcgr", self.sample + ".pcgr_acmg.grch38.maf")
            try:
                if self.ts_pcgr_maf != getime(pcgr_maf):
                    self.ts_pcgr_maf = getime(pcgr_maf)
                    self.pcgr_maf = pcgr_maf
            except FileNotFoundError:
                pass
            pcgr_snvs_indels = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble", self.sample, "pcgr", self.sample + ".pcgr_acmg.grch38.snvs_indels.tiers.tsv")
            try:
                if self.ts_pcgr_snvs_indels != getime(pcgr_snvs_indels):
                    self.ts_pcgr_snvs_indels = getime(pcgr_snvs_indels)
                    self.pcgr_snvs_indels = pcgr_snvs_indels
            except FileNotFoundError:
                pass
            pcgr_cna_segments = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble", self.sample, "pcgr", self.sample + ".pcgr_acmg.grch38.cna_segments.tsv.gz")
            try:
                if self.ts_pcgr_cna_segments != getime(pcgr_cna_segments):
                    self.ts_pcgr_cna_segments = getime(pcgr_cna_segments)
                    self.pcgr_cna_segments = pcgr_cna_segments
            except FileNotFoundError:
                pass

        if self.rna_true == "NA":
            self.rna_pcgr_report = "NA"
            self.rna_pcgr_maf = "NA"
            self.rna_pcgr_snvs_indels = "NA"
            self.rna_pcgr_cna_segments = "NA"
            self.ts_rna_pcgr_report = "NA"
            self.ts_rna_pcgr_maf = "NA"
            self.ts_rna_pcgr_snvs_indels = "NA"
            self.ts_rna_pcgr_cna_segments = "NA"
        else:
            rna_pcgr_report = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, "pcgr", self.rna + ".pcgr_acmg.grch38.flexdb.html")
            try:
                if self.ts_rna_pcgr_report != getime(rna_pcgr_report):
                    self.ts_rna_pcgr_report = getime(rna_pcgr_report)
                    self.rna_pcgr_report = rna_pcgr_report
            except FileNotFoundError:
                pass
            rna_pcgr_maf = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, "pcgr", self.rna + ".pcgr_acmg.grch38.maf")
            try:
                if self.ts_rna_pcgr_maf != getime(rna_pcgr_maf):
                    self.ts_rna_pcgr_maf = getime(rna_pcgr_maf)
                    self.rna_pcgr_maf = rna_pcgr_maf
            except FileNotFoundError:
                pass
            rna_pcgr_snvs_indels = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, "pcgr", self.rna + ".pcgr_acmg.grch38.snvs_indels.tiers.tsv")
            try:
                if self.ts_rna_pcgr_snvs_indels != getime(rna_pcgr_snvs_indels):
                    self.ts_rna_pcgr_snvs_indels = getime(rna_pcgr_snvs_indels)
                    self.rna_pcgr_snvs_indels = rna_pcgr_snvs_indels
            except FileNotFoundError:
                pass
            rna_pcgr_cna_segments = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, "pcgr", self.rna + ".pcgr_acmg.grch38.cna_segments.tsv.gz")
            try:
                if self.ts_rna_pcgr_cna_segments != getime(rna_pcgr_cna_segments):
                    self.ts_rna_pcgr_cna_segments = getime(rna_pcgr_cna_segments)
                    self.rna_pcgr_cna_segments = rna_pcgr_cna_segments
            except FileNotFoundError:
                pass

    def Gather_svariants(self):
        if self.dna_t_true == "NA":
            self.gridss = "NA"
            self.gripss_somatic = "NA"
            self.gripss_germline = "NA"
            self.purple_somatic = "NA"
            self.purple_germline = "NA"
            self.purple_circos = "NA"
            self.ts_gridss = "NA"
            self.ts_gripss_somatic = "NA"
            self.ts_gripss_germline = "NA"
            self.ts_purple_somatic = "NA"
            self.ts_purple_germline = "NA"
            self.ts_purple_circos = "NA"
        else:
            gridss = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/SVariants", self.sample, "gridss", self.dna_t + ".gridss.vcf.gz")
            try:
                if self.ts_gridss != getime(gridss):
                    self.ts_gridss = getime(gridss)
                    self.gridss = gridss
            except FileNotFoundError:
                pass
            gripss_somatic = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/SVariants", self.sample, "gridss", self.dna_t + ".gripss.filtered.somatic.vcf.gz")
            try:
                if self.ts_gripss_somatic != getime(gripss_somatic):
                    self.ts_gripss_somatic = getime(gripss_somatic)
                    self.gripss_somatic = gripss_somatic
            except FileNotFoundError:
                pass
            gripss_germline = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/SVariants", self.sample, "gridss", self.dna_n + ".gripss.filtered.germline.vcf.gz")
            try:
                if self.ts_gripss_germline != getime(gripss_germline):
                    self.ts_gripss_germline = getime(gripss_germline)
                    self.gripss_germline = gripss_germline
            except FileNotFoundError:
                pass
            purple_somatic = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/SVariants", self.sample, "purple", self.dna_t + ".driver.catalog.somatic.tsv")
            try:
                if self.ts_purple_somatic != getime(purple_somatic):
                    self.ts_purple_somatic = getime(purple_somatic)
                    self.purple_somatic = purple_somatic
            except FileNotFoundError:
                pass
            purple_germline = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/SVariants", self.sample, "purple", self.dna_t + ".driver.catalog.germline.tsv")
            try:
                if self.ts_purple_germline != getime(purple_germline):
                    self.ts_purple_germline = getime(purple_germline)
                    self.purple_germline = purple_germline
            except FileNotFoundError:
                pass
            purple_circos = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/SVariants", self.sample, "purple", "plot", self.dna_t + "circos.png")
            try:
                if self.ts_purple_circos != getime(purple_circos):
                    self.ts_purple_circos = getime(purple_circos)
                    self.purple_circos = purple_circos
            except FileNotFoundError:
                pass


    def Gather_RNA_other(self):
        if self.rna_true == "NA":
            self.annofuse = "NA"
            self.rna_abundance= "NA"
            self.big_wig_tracks_f = "NA"
            self.big_wig_tracks_r = "NA"
            self.ts_annofuse = "NA"
            self.ts_rna_abundance= "NA"
            self.ts_big_wig_tracks_f = "NA"
            self.ts_big_wig_tracks_r = "NA"
        else:
            annofuse_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/fusion", self.rna, "annoFuse", self.rna + ".putative_driver_fusions.tsv")
            try:
                if self.ts_annofuse != getime(annofuse_file):
                    self.ts_annofuse = getime(annofuse_file)
                    self.annofuse = annofuse_file
            except FileNotFoundError:
                pass

            # Change from stringtie to kallisto: DONE
            rna_abundance_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/kallisto", self.rna, "abundance_transcripts.tsv")
            try:
                if self.ts_rna_abundance != getime(rna_abundance_file):
                    self.ts_rna_abundance = getime(rna_abundance_file)
                    self.rna_abundance = rna_abundance_file
                    self.Gather_RNA_Light_ini()
            except FileNotFoundError:
                pass

            big_wig_tracks_f_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/tracks/bigWig", self.rna + ".forward.bw")
            try:
                if self.big_wig_tracks_f != getime(big_wig_tracks_f_file):
                    self.ts_big_wig_tracks_f = getime(big_wig_tracks_f_file)
                    self.big_wig_tracks_f = big_wig_tracks_f_file
            except FileNotFoundError:
                pass
            big_wig_tracks_r_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/tracks/bigWig", self.rna + ".reverse.bw")
            try:
                if self.big_wig_tracks_r != getime(big_wig_tracks_r_file):
                    self.ts_big_wig_tracks_r = getime(big_wig_tracks_r_file)
                    self.big_wig_tracks_r = big_wig_tracks_r_file
            except FileNotFoundError:
                pass


def getime(path):
    date = datetime.datetime.fromtimestamp(os.path.getmtime(path))
    return date.strftime("%Y/%m/%d")

def main():
    connection = create_connection("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db")
    patients = extract_sample_names(connection)
    print("Updating Database...")
    with progressbar.ProgressBar(max_value=len(patients), widgets=WIDGETS) as progress:
        for index, patient in enumerate(patients, 1):
            sample = Progress(connection, patient)
            update_timestamp_details(sample)
            update_fileloc_details(sample)
            sample.Gather_Run_Proc_BAM()
            sample.Gather_BAM_loc()
            sample.Gather_Final_BAMs()
            sample.Gather_VCFs()
            sample.Gather_PCGR()
            sample.Gather_svariants()
            sample.Gather_reports()
            sample.Gather_RNA_other()
            update_timestamp_details(sample)
            update_fileloc_details(sample)
            sample.Update_status()
            progress.update(index)
    print("Committing changes to Database...")
    connection.commit()
    connection.close()
    print("...Done.")


if __name__ == '__main__':
    main()
    #Update db with the objects
