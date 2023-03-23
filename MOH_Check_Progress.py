#!/usr/bin/env python3
import progressbar
import glob
import sys
import re
import os.path
import time
import datetime
import errno
from  DB_OPS import update_metrics_db,create_connection,extract_sample_details,extract_fileloc_details,extract_timestamp_details,update_timestamp_details,update_fileloc_details,extract_sample_names,update_status_db


WIDGETS = [' [', progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') - ', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']

class SampleData:
    def __init__(self, connection, sample):
        data = []
        self.conn = connection
        data = extract_sample_details(connection, sample)
        self.sample = data[0]
        self.sample_true = data[1]
        self.institution = data[2]
        self.cohort  = data[3]
        self.dna_n  = data[4]
        self.dna_n_true = data[5]
        self.dna_t = data[6]
        self.dna_t_true = data[7]
        self.rna = data[8]
        self.rna_true = data[9]
    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

class Progress(SampleData):
    def grab_db_values(self):
        data = []
        data = extract_fileloc_details(self.conn, self.sample)
        self.run_proc_bam_dna_t = data[1]
        self.run_proc_bam_dna_n = data[2]
        self.beluga_bam_dna_t = data[3]
        self.beluga_bam_dna_n = data[4]
        self.dna_vcf_g = data[5]
        self.dna_vcf_s = data[6]
        self.mutect2_somatic_vcf = data[7]
        self.mutect2_germline_vcf = data[8]
        self.strelka2_germline_vcf = data[9]
        self.strelka2_somatic_vcf = data[10]
        self.vardict_germline_vcf = data[11]
        self.vardict_somatic_vcf = data[12]
        self.varscan2_germline_vcf = data[13]
        self.varscan2_somatic_vcf = data[14]
        self.cnvkit_vcf = data[15]
        self.final_vcf = data[16]
        self.final_dna_bam_t = data[17]
        self.final_dna_bam_n = data[18]
        self.dna_multiqc = data[19]
        self.pcgr_report = data[20]
        self.pcgr_maf = data[21]
        self.pcgr_snvs_indels = data[22]
        self.pcgr_cna_segments = data[23]
        self.tp_ini = data[24]
        self.run_proc_fastq_1_rna = data[25]
        self.run_proc_fastq_2_rna = data[26]
        self.beluga_fastq_1_rna = data[27]
        self.beluga_fastq_2_rna = data[28]
        self.rna_vcf = data[29]
        # self.final_rna_bam_expression = data[26]
        self.final_rna_bam_variants = data[30]
        self.rna_multiqc = data[31]
        self.annofuse = data[32]
        self.gridss = data[33]
        self.rna_abundance = data[34]
        self.big_wig_tracks_f = data[35]
        self.big_wig_tracks_r = data[36]
        self.rna_abundance_ini = data[37]
        self.rna_variants_ini = data[38]
        #get the old timestamps
        data = []
        data = extract_timestamp_details(self.conn, self.sample)
        self.ts_run_proc_bam_dna_t = data[1]
        self.ts_run_proc_bam_dna_n = data[2]
        self.ts_beluga_bam_dna_t = data[3]
        self.ts_beluga_bam_dna_n = data[4]
        self.ts_dna_vcf_g = data[5]
        self.ts_dna_vcf_s = data[6]
        self.ts_mutect2_somatic_vcf = data[7]
        self.ts_mutect2_germline_vcf = data[8]
        self.ts_strelka2_germline_vcf = data[9]
        self.ts_strelka2_somatic_vcf = data[10]
        self.ts_vardict_germline_vcf = data[11]
        self.ts_vardict_somatic_vcf = data[12]
        self.ts_varscan2_germline_vcf = data[13]
        self.ts_varscan2_somatic_vcf = data[14]
        self.ts_cnvkit_vcf = data[15]
        self.ts_final_vcf = data[16]
        self.ts_final_dna_bam_t = data[17]
        self.ts_final_dna_bam_n = data[18]
        self.ts_dna_multiqc = data[19]
        self.ts_pcgr_report = data[20]
        self.ts_pcgr_maf = data[21]
        self.ts_pcgr_snvs_indels = data[22]
        self.ts_pcgr_cna_segments = data[23]
        self.ts_tp_ini = data[24]
        self.ts_run_proc_fastq_1_rna = data[25]
        self.ts_run_proc_fastq_2_rna = data[26]
        self.ts_beluga_fastq_1_rna = data[27]
        self.ts_beluga_fastq_2_rna = data[28]
        self.ts_rna_vcf = data[29]
        # self.ts_final_rna_bam_expression = data[26]
        self.ts_final_rna_bam_variants = data[30]
        self.ts_rna_multiqc = data[31]
        self.ts_annofuse = data[32]
        self.ts_gridss = data[33]
        self.ts_rna_abundance = data[34]
        self.ts_big_wig_tracks_f = data[35]
        self.ts_big_wig_tracks_r = data[36]
        self.ts_rna_abundance_ini = data[37]
        self.ts_rna_variants_ini = data[38]

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
        dna_n_transfered = "NA"
        if self.ts_run_proc_bam_dna_n != "NA":
            dna_n_transfered = "Complete"
        dna_t_transfered = "NA"
        if self.ts_run_proc_bam_dna_t != "NA":
            dna_t_transfered = "Complete"
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
        if dna_n_transfered == dna_t_transfered == dna_alignment == dna_variant_call == dna_report == "Complete":
            dna_pipeline_execution = "Complete"
        # RNA
        rna_transferred = "NA"
        if self.beluga_fastq_1_rna != "NA":
            rna_transferred = "Complete"
        # For abundance
        rna_pipeline_light_execution = "NA"
        if self.rna_abundance != "NA":
            rna_pipeline_light_execution = "Complete"
        rna_pseudoalignment = "NA"
        # if self.final_rna_bam_expression != "NA":
        #     rna_pseudoalignment = "Complete"
        rna_variant_call = "NA"
        if self.final_rna_bam_variants != "NA":
            rna_variant_call = "Complete"
        rna_report = "NA"
        if self.annofuse != "NA":
            rna_report = "Complete"
        rna_pipeline_execution = "NA"
        if rna_transferred == rna_pipeline_light_execution == rna_pseudoalignment == rna_variant_call == rna_report == "Complete":
            rna_pipeline_execution = "Complete"
        # DELIVERY
        dna_delivered, rna_light_delivered, rna_delivered = self.check_delivery()
        update_status_db(
            self.conn,
            self.sample,
            dna_n_transfered,
            dna_t_transfered,
            dna_alignment,
            dna_variant_call,
            dna_report,
            dna_pipeline_execution,
            dna_delivered,
            rna_transferred,
            rna_pipeline_light_execution,
            rna_light_delivered,
            rna_pseudoalignment,
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
                                self.run_proc_bam_dna_n = fields[0]
                                self.ts_run_proc_bam_dna_n  = getime(filename)
                                dna_n = True
                            elif self.dna_t in line:
                                fields = line.split(",")
                                self.run_proc_bam_dna_t = fields[0]
                                self.ts_run_proc_bam_dna_t = getime(filename)
                                dna_t = True

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
                                    self.run_proc_fastq_1_rna = fields[0]
                                    self.ts_run_proc_fastq_1_rna = getime(filename)
                                    fastq1 = True
                                elif re.search(r"R2.*.fastq", os.path.basename(fields[0])):
                                    self.run_proc_fastq_2_rna = fields[0]
                                    self.ts_run_proc_fastq_2_rna = getime(filename)
                                    fastq2 = True

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
            # else:
            #     raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"{dna_n_file} nor {dna_n_file_old}")
            # bam dna_t in raw_reads
            dna_t_file_rr = os.path.join(loc1, self.dna_t, os.path.basename(self.run_proc_bam_dna_t))
            dna_t_file_old_rr = os.path.join(loc1, self.dna_t, self.dna_t + ".bam")
            if os.path.exists(dna_t_file_rr):
                self.beluga_bam_dna_t = dna_t_file_rr
                self.ts_beluga_bam_dna_t = getime(dna_t_file_rr)
            elif os.path.exists(dna_t_file_old_rr):
                self.beluga_bam_dna_t = dna_t_file_old_rr
                self.ts_beluga_bam_dna_t = getime(dna_t_file_old_rr)
            # else:
            #     raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"{dna_t_file} nor {dna_t_file_old}")
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
            # else:
            #     raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"{rna_fq1} nor {rna_fq1_old}")
            # fastq2 in raw_reads
            rna_fq2_rr = os.path.join(loc1, self.rna, os.path.basename(self.run_proc_fastq_2_rna))
            rna_fq2_old_rr = os.path.join(loc1, self.rna, self.rna + "_R2.fastq.gz")
            if os.path.exists(rna_fq2_rr):
                self.beluga_fastq_2_rna = rna_fq2_rr
                self.ts_beluga_fastq_2_rna = getime(rna_fq2_rr)
            elif os.path.exists(rna_fq2_old_rr):
                self.beluga_fastq_2_rna = rna_fq2_old_rr
                self.ts_beluga_fastq_2_rna = getime(rna_fq2_old_rr)
            # else:
            #     raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"{rna_fq2} nor {rna_fq2_old}")
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
            # for filename in glob.glob(path):
            #     if self.ts_final_dna_bam_n != getime(filename):
            #         self.final_dna_bam_n = filename
            #         self.ts_final_dna_bam_n = getime(filename)
            #         self.Gather_DNA_ini()
            final_dna_bam_t_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.dna_t, self.dna_t + ".sorted.dup.recal.bam")
            try:
                if self.ts_final_dna_bam_t != getime(final_dna_bam_t_file):
                    self.final_dna_bam_t = final_dna_bam_t_file
                    self.ts_final_dna_bam_t = getime(final_dna_bam_t_file)
                    self.Gather_DNA_ini()
            except FileNotFoundError:
                pass
            # path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment/" + self.dna_t + "/*T.sorted.dup.recal.bam"
            # for filename in glob.glob(path):
            #     if self.ts_final_dna_bam_t != getime(filename):
            #         self.final_dna_bam_t = filename
            #         self.ts_final_dna_bam_t = getime(filename)
            #         self.Gather_DNA_ini()
        if self.rna_true == "NA":
            # self.final_rna_bam_expression = "NA"
            self.final_rna_bam_variants = "NA"
            # self.ts_final_rna_bam_expression = "NA"
            self.ts_final_rna_bam_variants = "NA"
        else:
            final_rna_bam_variants_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, self.rna + ".sorted.mdup.cram")
            try:
                if self.ts_final_rna_bam_variants != getime(final_rna_bam_variants_file):
                    self.ts_final_rna_bam_variants = getime(final_rna_bam_variants_file)
                    self.final_rna_bam_expression = final_rna_bam_variants_file
            except FileNotFoundError:
                pass
            # final_rna_bam_expression_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, self.rna + ".sorted.mdup.split.realigned.recal.bam")
            # try:
            #     if self.ts_final_rna_bam_expression != getime(final_rna_bam_expression_file):
            #         self.ts_final_rna_bam_expression = getime(final_rna_bam_expression_file)
            #         self.final_rna_bam_variants = final_rna_bam_expression_file
            # except FileNotFoundError:
            #     pass
            # path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment/" + self.rna + "*/*"
            # for filename in glob.glob(path):
            #     if ".sorted.mdup.cram" in filename:
            #         if self.ts_final_rna_bam_variants != getime(filename):
            #             self.ts_final_rna_bam_variants = getime(filename)
            #             self.final_rna_bam_expression = filename
            #     elif ".sorted.mdup.split.realigned.recal.bam" in filename:
            #         if ".bai" in filename or ".sbi" in filename:
            #             continue
            #         if self.ts_final_rna_bam_expression != getime(filename):
            #             self.ts_final_rna_bam_expression = getime(filename)
            #             self.final_rna_bam_variants = filename

    def Gather_VCFs(self):
        if self.dna_n_true == "NA":
            self.dna_vcf_g = "NA"
            self.dna_vcf_s = "NA"
            self.mutect2_somatic_vcf = "NA"
            self.mutect2_germline_vcf = "NA"
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
            self.ts_mutect2_germline_vcf = "NA"
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
            # path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedvariants/ensemble/" + self.sample + "*/*.vt.annot.vcf.gz"
            # for filename in glob.glob(path):
            #     if "ensemble.germline.vt.annot.vcf.gz" in filename:
            #         if self.ts_dna_vcf_g != getime(filename):
            #             self.ts_dna_vcf_g = getime(filename)
            #             self.dna_vcf_g = filename
            #     elif "ensemble.somatic.vt.annot.vcf.gz" in filename:
            #         if self.ts_dna_vcf_s != getime(filename):
            #             self.ts_dna_vcf_s = getime(filename)
            #             self.dna_vcf_s = filename
            mutect2_somatic_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants", self.sample, self.sample + ".mutect2.somatic.vt.vcf.gz")
            try:
                if self.ts_mutect2_somatic_vcf != getime(mutect2_somatic_vcf_file):
                    self.ts_mutect2_somatic_vcf = getime(mutect2_somatic_vcf_file)
                    self.mutect2_somatic_vcf = mutect2_somatic_vcf_file
            except FileNotFoundError:
                pass
            mutect2_germline_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants", self.sample, self.sample + ".mutect2.vcf.gz")
            try:
                if self.ts_mutect2_germline_vcf != getime(mutect2_germline_vcf_file):
                    self.ts_mutect2_germline_vcf = getime(mutect2_germline_vcf_file)
                    self.mutect2_germline_vcf = mutect2_germline_vcf_file
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
            # path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedvariants/" + self.sample + "*/*.vcf.gz"
            # for filename in glob.glob(path):
            #     if ".mutect2.somatic.vt.vcf.gz" in filename:
            #         if self.ts_mutect2_somatic_vcf != getime(filename):
            #             self.ts_mutect2_somatic_vcf = getime(filename)
            #             self.mutect2_somatic_vcf = filename
            #     if ".mutect2.vcf.gz" in filename:
            #         if self.ts_mutect2_germline_vcf != getime(filename):
            #             self.ts_mutect2_germline_vcf = getime(filename)
            #             self.mutect2_germline_vcf = filename
            #     if ".strelka2.somatic.purple.vcf.gz" in filename:
            #         if self.ts_strelka2_somatic_vcf != getime(filename):
            #             self.ts_strelka2_somatic_vcf = getime(filename)
            #             self.strelka2_somatic_vcf = filename
            #     if ".strelka2.germline.vt.vcf.gz" in filename:
            #         if self.ts_strelka2_germline_vcf != getime(filename):
            #             self.ts_strelka2_germline_vcf = getime(filename)
            #             self.strelka2_germline_vcf = filename
            #     if ".vardict.germline.vt.vcf.gz" in filename:
            #         if self.ts_vardict_germline_vcf != getime(filename):
            #             self.ts_vardict_germline_vcf = getime(filename)
            #             self.vardict_germline_vcf= filename
            #     if ".vardict.somatic.vt.vcf.gz" in filename:
            #         if self.ts_vardict_somatic_vcf != getime(filename):
            #             self.ts_vardict_somatic_vcf = getime(filename)
            #             self.vardict_somatic_vcf = filename
            #     if ".varscan2.germline.vt.vcf.gz" in filename:
            #         if self.ts_varscan2_germline_vcf != getime(filename):
            #             self.ts_varscan2_germline_vcf = getime(filename)
            #             self.varscan2_germline_vcf = filename
            #     if ".varscan2.somatic.vt.vcf.gz" in filename:
            #         if self.ts_varscan2_somatic_vcf != getime(filename):
            #             self.ts_varscan2_somatic_vcf = getime(filename)
            #             self.varscan2_somatic_vcf = filename
        if self.rna_true == "NA":
            self.rna_vcf = "NA"
            self.ts_rna_vcf= "NA"
        else:
            rna_vcf_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment", self.rna, self.rna + ".hc.vcf.gz")
            try:
                if self.ts_rna_vcf != getime(rna_vcf_file):
                    self.rna_vcf = rna_vcf_file
                    self.ts_rna_vcf = getime(rna_vcf_file)
                    self.Gather_RNA_Variants_ini()
            except FileNotFoundError:
                pass
            # path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment/" + self.rna + "/*.hc.vcf.gz"
            # for filename in glob.glob(path):
            #     if self.ts_rna_vcf != getime(filename):
            #         self.rna_vcf = filename
            #         self.ts_rna_vcf = getime(filename)
            #         self.Gather_RNA_ini()


    def Gather_reports(self):
        if self.dna_n_true == "NA":
            self.dna_multiqc = "NA"
            # self.pcgr = "NA"
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
            # path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna/" + self.sample + "*multiqc.html"
            # for filename in glob.glob(path):
            #     if self.ts_dna_multiqc != getime(filename):
            #         self.ts_dna_multiqc = getime(filename)
            #         self.dna_multiqc = filename

            # "[sample].pcgr_acmg.grch38.flexdb.html"
            # "[sample].pcgr_acmg.grch38.maf"
            # "[sample].pcgr_acmg.grch38.snvs_indels.tiers.tsv"
            # "[sample].pcgr_acmg.grch38.cna_segments.tsv.gz"
            # "raw cnv ici : /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/SVariants/[sample]/[sample].cnvkit.vcf.gz"
            # pcgr_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble", self.sample, "pcgr", self.sample + ".pcgr_acmg.grch38.flexdb.html")
            # try:
            #     if self.ts_pcgr != getime(pcgr_file):
            #         self.ts_pcgr = getime(pcgr_file)
            #         self.pcgr = pcgr_file
            # except FileNotFoundError:
            #     pass
            # path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedvariants/ensemble/" + self.sample + "*/pcgr*/*.flexdb.html"
            # for filename in glob.glob(path):
            #     if self.ts_pcgr != getime(filename):
            #         self.ts_pcgr = getime(filename)
            #         self.pcgr = filename

        if self.rna_true == "NA":
            self.rna_multiqc = "NA"
            self.ts_rna_multiqc = "NA"
        else:
            #####################Not Implemented########################################
            self.rna_multiqc = "NA"
            self.ts_rna_multiqc = "NA"
            ############################################################################

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

    def Gather_RNA_other(self):
        if self.rna_true == "NA":
            self.annofuse = "NA"
            self.gridss= "NA"
            self.rna_abundance= "NA"
            self.big_wig_tracks_f = "NA"
            self.big_wig_tracks_r = "NA"
            self.ts_annofuse = "NA"
            self.ts_gridss= "NA"
            self.ts_rna_abundance= "NA"
            self.ts_big_wig_tracks_f = "NA"
            self.ts_big_wig_tracks_r = "NA"
        else:
            #####################Not Implemented########################################
            self.gridss= "NA"
            self.ts_gridss= "NA"
            ############################################################################
            annofuse_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/fusion", self.rna, "annoFuse", self.rna + ".putative_driver_fusions.tsv")
            try:
                if self.ts_annofuse != getime(annofuse_file):
                    self.ts_annofuse = getime(annofuse_file)
                    self.annofuse = annofuse_file
            except FileNotFoundError:
                pass
            # path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/fusion/" + self.rna + "/annoFuse/*.putative_driver_fusions.tsv"
            # for filename in glob.glob(path):
            #     if self.ts_annofuse != getime(filename):
            #         self.ts_annofuse = getime(filename)
            #         self.annofuse = filename

            # Change from stringtie to kallisto: DONE
            rna_abundance_file = os.path.join("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/kallisto", self.rna, "abundance_transcripts.tsv")
            try:
                if self.ts_rna_abundance != getime(rna_abundance_file):
                    self.ts_rna_abundance = getime(rna_abundance_file)
                    self.rna_abundance = rna_abundance_file
                    self.Gather_RNA_Light_ini()
            except FileNotFoundError:
                pass
            # for filename in glob.glob(path):
            #     if self.ts_rna_abundance != getime(filename):
            #         self.ts_rna_abundance = getime(filename)
            #         self.rna_abundance = filename

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
            # path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/tracks/bigWig/" + self.rna + "*"
            # for filename in glob.glob(path):
            #     if "forward" in filename:
            #         if self.big_wig_tracks_f != getime(filename):
            #             self.ts_big_wig_tracks_f = getime(filename)
            #             self.big_wig_tracks_f = filename
            #     if "reverse" in filename:
            #         if self.big_wig_tracks_r != getime(filename):
            #             self.ts_big_wig_tracks_r = getime(filename)
            #             self.big_wig_tracks_r = filename


def getime(path):
    date = datetime.datetime.fromtimestamp(os.path.getmtime(path))
    return date.strftime("%Y/%m/%d")

def main():
    connection = create_connection("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db")
    # ALL_Samples = extract_sample_names(connection)
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
