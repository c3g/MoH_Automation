#!/usr/bin/env python3
import glob
import sys
import re
import os.path
import time
import datetime
from  DB_OPS import update_metrics_db,create_connection,extract_sample_details,extract_fileloc_details,extract_timestamp_details,update_timestamp_details,update_fileloc_details,extract_sample_names,update_status_db


class SampleData:
    def __init__(self, connection, Sample):
        data = []
        self.conn = connection
        data = extract_sample_details(connection,Sample) 
        self.Sample = data[0]
        self.Sample_True = data[1]
        self.Institution = data[2]
        self.Cohort  = data[3]
        self.DNA_N  = data[4]
        self.DNA_N_True = data[5]
        self.DNA_T = data[6]
        self.DNA_T_True = data[7]
        self.RNA = data[8]
        self.RNA_True = data[9]
    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

class Progress(SampleData):
    def grab_db_values(self):
        data = []
        data = extract_fileloc_details(self.conn,self.Sample) 
        self.Run_Proc_BAM_DNA_T = data[1]
        self.Run_Proc_BAM_DNA_N = data[2]
        self.Beluga_BAM_DNA_T = data[3]
        self.Beluga_BAM_DNA_N = data[4]
        self.DNA_VCF_G = data[5]
        self.DNA_VCF_S = data[6]
        self.Mutect2_Somatic_vcf = data[7]
        self.Mutect2_Germline_vcf = data[8]
        self.strelka2_Germline_vcf = data[9]
        self.strelka2_Somatic_vcf = data[10]
        self.vardict_Germline_vcf = data[11]
        self.vardict_Somatic_vcf = data[12]
        self.varscan2_Germline_vcf = data[13]
        self.varscan2_Somatic_vcf = data[14]
        self.Final_VCF = data[15]
        self.Final_DNA_BAM_T = data[16]
        self.Final_DNA_BAM_N = data[17]
        self.DNA_MultiQC = data[18]
        self.PCGR = data[19]
        self.TP_ini = data[20]
        self.Run_Proc_fastq_1_RNA = data[21]
        self.Run_Proc_fastq_2_RNA = data[22]
        self.Beluga_fastq_1_RNA = data[23]
        self.Beluga_fastq_2_RNA = data[24]
        self.RNA_VCF = data[25]
        self.Final_RNA_BAM_expression = data[26]
        self.Final_RNA_BAM_variants = data[27]
        self.RNA_MultiQC = data[28]
        self.AnnoFuse = data[29]
        self.GRIDSS = data[30]
        self.RNA_Abundance = data[31]
        self.big_wig_tracks_F = data[32]
        self.big_wig_tracks_R = data[33]
        self.RNA_Abundance_ini = data[34]
        self.RNA_Variants_ini = data[35]
        #get the old timestamps
        data = []
        data = extract_timestamp_details(self.conn,self.Sample) 
        self.TS_Run_Proc_BAM_DNA_T = data[1]
        self.TS_Run_Proc_BAM_DNA_N = data[2]
        self.TS_Beluga_BAM_DNA_T = data[3]
        self.TS_Beluga_BAM_DNA_N = data[4]
        self.TS_DNA_VCF_G = data[5]
        self.TS_DNA_VCF_S = data[6]
        self.TS_Mutect2_Somatic_vcf = data[7]
        self.TS_Mutect2_Germline_vcf = data[8]
        self.TS_strelka2_Germline_vcf = data[9]
        self.TS_strelka2_Somatic_vcf = data[10]
        self.TS_vardict_Germline_vcf = data[11]
        self.TS_vardict_Somatic_vcf = data[12]
        self.TS_varscan2_Germline_vcf = data[13]
        self.TS_varscan2_Somatic_vcf = data[14]
        self.TS_Final_VCF = data[15]
        self.TS_Final_DNA_BAM_T = data[16]
        self.TS_Final_DNA_BAM_N = data[17]
        self.TS_DNA_MultiQC = data[18]
        self.TS_PCGR = data[19]
        self.TS_TP_ini = data[20]
        self.TS_Run_Proc_fastq_1_RNA = data[21]
        self.TS_Run_Proc_fastq_2_RNA = data[22]
        self.TS_Beluga_fastq_1_RNA = data[23]
        self.TS_Beluga_fastq_2_RNA = data[24]
        self.TS_RNA_VCF = data[25]
        self.TS_Final_RNA_BAM_expression = data[26]
        self.TS_Final_RNA_BAM_variants = data[27]
        self.TS_RNA_MultiQC = data[28]
        self.TS_AnnoFuse = data[29]
        self.TS_GRIDSS = data[30]
        self.TS_RNA_Abundance = data[31]
        self.TS_big_wig_tracks_F = data[32]
        self.TS_big_wig_tracks_R = data[33]
        self.TS_RNA_Abundance_ini = data[34]
        self.TS_RNA_Variants_ini = data[35]

    def __init__(self, connection, Sample):
        super().__init__(connection, Sample)
        self.grab_db_values()

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

    def Gather_DNA_ini(self):
#################ONCE DONE JUST STORE THE WHOLE FLIPPEN INI HERE ########################################
        filename = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/TumorPair.config.trace.ini'
        self.TP_ini = filename
        self.TS_TP_ini = getime(filename)
#########################################################################################################
    def Gather_RNA_ini(self):
#################ONCE DONE JUST STORE THE WHOLE FLIPPEN INI HERE########################################
        filename = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/RnaSeq.config.trace.ini'
        self.RNA_Abundance_ini = filename
        self.TS_RNA_Abundance_ini = getime(filename)
        filename = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/RnaSeq.config.trace.ini'
        self.RNA_Variants_ini = filename
        self.TS_RNA_Variants_ini = getime(filename)
#########################################################################################################

    def Update_status(self):
        DNA_N_Transfered = "NA"
        if self.TS_Run_Proc_BAM_DNA_N != "NA":
            DNA_N_Transfered = "Complete"
        DNA_T_Transfered = "NA"
        if self.TS_Run_Proc_BAM_DNA_T != "NA":
            DNA_T_Transfered = "Complete"
        Alignment= "NA"    
        if self.Final_DNA_BAM_T != "NA":
            Alignment = "Complete"
        Variants = "NA"    
        if self.DNA_VCF_G != "NA" and self.DNA_VCF_S != "NA" :
            Variants = "Complete"
        Reports = "NA"    
        if self.DNA_MultiQC != "NA" and self.PCGR != "NA":
            Reports = "Complete"
        Tumour_Pair_Complete = "NA"
        if Reports != "NA" and Alignment != "NA" and Variants != "NA":
            Tumour_Pair_Complete = "Complete"
        RNA_Transfrered = "NA"
        if self.Beluga_fastq_1_RNA != "NA":
            RNA_Transfrered = "Complete"
        RNA_Alignment_expression = "NA"
        if self.Final_RNA_BAM_expression != "NA":
            RNA_Alignment_expression = "Complete"
        RNA_Alignment_Variant = "NA"
        if self.Final_RNA_BAM_variants != "NA":
            RNA_Alignment_Variant = "Complete"
        RNA_Reports = "NA"
        if self.AnnoFuse != "NA":
            RNA_Reports = "Complete"
        RNA_Complete = "NA"
        if self.big_wig_tracks_R != "NA" and RNA_Alignment_expression != "NA" and RNA_Alignment_Variant != "NA" and RNA_Reports != "NA":
            RNA_Complete = "Complete"
        update_status_db(self.conn,self.Sample,DNA_N_Transfered,DNA_T_Transfered,Alignment,Variants,Reports,Tumour_Pair_Complete,RNA_Transfrered,RNA_Alignment_expression,RNA_Alignment_Variant,RNA_Reports,RNA_Complete)


    def Gather_Run_Proc_BAM(self):
        if self.DNA_N_True == "NA":
            self.Run_Proc_BAM_DNA_T = "NA"
            self.Run_Proc_BAM_DNA_N = "NA"
            self.TS_Run_Proc_BAM_DNA_T = "NA"
            self.TS_Run_Proc_BAM_DNA_N = "NA"
        else:
            path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*'
            if self.Run_Proc_BAM_DNA_T == "NA" or self.Run_Proc_BAM_DNA_N == "NA":
                for filename in glob.glob(path):
                    with open(filename, 'r') as f:
                        for line in f:
                            if self.DNA_N in line:
                               fields = line.split(",")
                               self.Run_Proc_BAM_DNA_N = fields[0] 
                               self.TS_Run_Proc_BAM_DNA_N  = getime(filename)
                            elif self.DNA_T in line:
                               fields = line.split(",")
                               self.Run_Proc_BAM_DNA_T = fields[0]
                               self.TS_Run_Proc_BAM_DNA_T = getime(filename)

        if self.RNA_True == "NA":
            self.Run_Proc_fastq_1_RNA = "NA"
            self.Run_Proc_fastq_2_RNA = "NA"
            self.TS_Run_Proc_fastq_1_RNA = "NA"
            self.TS_Run_Proc_fastq_2_RNA = "NA"
        else:
            path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/log_files/transfer/*'
            if self.Run_Proc_fastq_1_RNA == "NA" or self.Run_Proc_fastq_2_RNA == "NA":
                for filename in glob.glob(path):
                    with open(filename, 'r') as f:
                        for line in f:
                            if self.RNA in line:
                                if "R1.fastq" in line:
                                    fields = line.split(",")
                                    self.Run_Proc_fastq_1_RNA = fields[0] 
                                    self.TS_Run_Proc_fastq_1_RNA = getime(filename)
                                elif "R2.fastq" in line:
                                    fields = line.split(",")
                                    self.Run_Proc_fastq_2_RNA = fields[0] 
                                    self.TS_Run_Proc_fastq_2_RNA = getime(filename)

    def Gather_BAM_loc(self):
        loc1 = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/raw_reads"
        loc2 = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/raw_reads"
        if self.DNA_N_True == "NA":
            self.Beluga_BAM_DNA_T = "NA"
            self. Beluga_BAM_DNA_N= "NA"
            self.TS_Beluga_BAM_DNA_T = "NA"
            self.TS_Beluga_BAM_DNA_N = "NA"
        else:
            path = loc1 + "/*" + self.DNA_N + "*/*.bam"
            for filename in glob.glob(path):
               self.Beluga_BAM_DNA_N = filename
               self.TS_Beluga_BAM_DNA_N = getime(filename)
            path = loc1 + "/*" + self.DNA_T + "*/*.bam"
            for filename in glob.glob(path):
               self.Beluga_BAM_DNA_T = filename
               self.TS_Beluga_BAM_DNA_T = getime(filename)
            path = loc2 + "/*" + self.DNA_N + "*/*.bam"
            for filename in glob.glob(path):
               self.Beluga_BAM_DNA_N = filename
               self.TS_Beluga_BAM_DNA_N = getime(filename)
            path = loc2 + "/*" + self.DNA_T + "*/*.bam"
            for filename in glob.glob(path):
               self.Beluga_BAM_DNA_T = filename
               self.TS_Beluga_BAM_DNA_T = getime(filename)
        if self.RNA_True == "NA":
            self.Beluga_fastq_1_RNA = "NA"
            self.Beluga_fastq_2_RNA = "NA"
            self.TS_Beluga_fastq_1_RNA = "NA"
            self.TS_Beluga_fastq_2_RNA = "NA"
        else:
            path = loc1 + "/*" + self.RNA+ "*/*.fastq.gz"
            for filename in glob.glob(path):
                if "R1.fastq.gz" in filename:
                    self.Beluga_fastq_1_RNA = filename
                    self.TS_Beluga_fastq_1_RNA = getime(filename)
                elif "R2.fastq.gz" in filename:
                    self.Beluga_fastq_2_RNA = filename
                    self.TS_Beluga_fastq_2_RNA = getime(filename)
            path = loc2 + "/*" + self.RNA+ "*/*.fastq.gz"
            for filename in glob.glob(path):
                if "R1.fastq.gz" in filename:
                    self.Beluga_fastq_1_RNA = filename
                    self.TS_Beluga_fastq_1_RNA = getime(filename)
                elif "R2.fastq.gz" in filename:
                    self.Beluga_fastq_2_RNA = filename
                    self.TS_Beluga_fastq_2_RNA = getime(filename)

    def Gather_Final_BAMs(self):
        if self.DNA_N_True == "NA":
            self.Final_DNA_BAM_T = "NA"
            self.Final_DNA_BAM_N = "NA"
            self.TS_Final_DNA_BAM_T= "NA"
            self.TS_Final_DNA_BAM_N= "NA"
        else:
            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment/" + self.DNA_N + "/*N.sorted.dup.recal.bam"
            for filename in glob.glob(path):
                if self.TS_Final_DNA_BAM_N != getime(filename):
                    self.Final_DNA_BAM_N = filename
                    self.TS_Final_DNA_BAM_N = getime(filename)
                    self.Gather_DNA_ini()
            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment/" + self.DNA_T + "/*T.sorted.dup.recal.bam"
            for filename in glob.glob(path):
                if self.TS_Final_DNA_BAM_T != getime(filename):
                    self.Final_DNA_BAM_T = filename
                    self.TS_Final_DNA_BAM_T = getime(filename)
                    self.Gather_DNA_ini()
        if self.RNA_True == "NA":
            self.Final_RNA_BAM_expression = "NA"
            self.Final_RNA_BAM_variants = "NA"
            self.TS_Final_RNA_BAM_expression = "NA"
            self.TS_Final_RNA_BAM_variants = "NA"
        else:
            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment/" + self.RNA + "*/*"
            for filename in glob.glob(path):
                if ".sorted.mdup.cram" in filename:
                    if self.TS_Final_RNA_BAM_variants != getime(filename):
                        self.TS_Final_RNA_BAM_variants = getime(filename)
                        self.Final_RNA_BAM_expression = filename
                elif ".sorted.mdup.split.realigned.recal.bam" in filename:
                    if ".bai" in filename or ".sbi" in filename:
                        continue
                    if self.TS_Final_RNA_BAM_expression != getime(filename):
                        self.TS_Final_RNA_BAM_expression = getime(filename)
                        self.Final_RNA_BAM_variants = filename

    def Gather_VCFs(self):
        if self.DNA_N_True == "NA":
            self.DNA_VCF_G = "NA"
            self.DNA_VCF_S = "NA"
            self.Mutect2_Somatic_vcf = "NA"
            self.Mutect2_Germline_vcf = "NA"
            self.strelka2_Germline_vcf = "NA"
            self.strelka2_Somatic_vcf = "NA"
            self.vardict_Germline_vcf = "NA"
            self.vardict_Somatic_vcf = "NA"
            self.varscan2_Germline_vcf = "NA"
            self.varscan2_Somatic_vcf = "NA"
            self.Final_VCF = "NA"
            self.TS_DNA_VCF_G = "NA"
            self.TS_DNA_VCF_S = "NA"
            self.TS_Mutect2_Somatic_vcf = "NA"
            self.TS_Mutect2_Germline_vcf = "NA"
            self.TS_strelka2_Germline_vcf = "NA"
            self.TS_strelka2_Somatic_vcf = "NA"
            self.TS_vardict_Germline_vcf = "NA"
            self.TS_vardict_Somatic_vcf = "NA"
            self.TS_varscan2_Germline_vcf = "NA"
            self.TS_varscan2_Somatic_vcf = "NA"
            self.TS_Final_VCF = "NA"
        else:
            #####################Not Implemented########################################
            self.TS_Final_VCF = "NA"
            self.Final_VCF = "NA"
            ############################################################################
            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble/" + self.Sample + "*/*.vt.annot.vcf.gz"
            for filename in glob.glob(path):
                if "ensemble.germline.vt.annot.vcf.gz" in filename:
                    if self.TS_DNA_VCF_G != getime(filename):
                        self.TS_DNA_VCF_G = getime(filename)
                        self.DNA_VCF_G = filename
                elif "ensemble.somatic.vt.annot.vcf.gz" in filename:
                    if self.TS_DNA_VCF_S != getime(filename):
                        self.TS_DNA_VCF_S = getime(filename)
                        self.DNA_VCF_S = filename
            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/" + self.Sample + "*/*.vcf.gz"
            for filename in glob.glob(path):
                if ".mutect2.somatic.vt.vcf.gz" in filename:
                    if self.TS_Mutect2_Somatic_vcf != getime(filename):
                        self.TS_Mutect2_Somatic_vcf = getime(filename)
                        self.Mutect2_Somatic_vcf = filename
                if ".mutect2.vcf.gz" in filename:
                    if self.TS_Mutect2_Germline_vcf != getime(filename):
                        self.TS_Mutect2_Germline_vcf = getime(filename)
                        self.Mutect2_Germline_vcf = filename
                if ".strelka2.somatic.purple.vcf.gz" in filename:
                    if self.TS_strelka2_Somatic_vcf != getime(filename):
                        self.TS_strelka2_Somatic_vcf = getime(filename)
                        self.strelka2_Somatic_vcf = filename
                if ".strelka2.germline.vt.vcf.gz" in filename:
                    if self.TS_strelka2_Germline_vcf != getime(filename):
                        self.TS_strelka2_Germline_vcf = getime(filename)
                        self.strelka2_Germline_vcf = filename
                if ".vardict.germline.vt.vcf.gz" in filename:
                    if self.TS_vardict_Germline_vcf != getime(filename):
                        self.TS_vardict_Germline_vcf = getime(filename)
                        self.vardict_Germline_vcf= filename
                if ".vardict.somatic.vt.vcf.gz" in filename:
                    if self.TS_vardict_Somatic_vcf != getime(filename):
                        self.TS_vardict_Somatic_vcf = getime(filename)
                        self.vardict_Somatic_vcf = filename
                if ".varscan2.germline.vt.vcf.gz" in filename:
                    if self.TS_varscan2_Germline_vcf != getime(filename):
                        self.TS_varscan2_Germline_vcf = getime(filename)
                        self.varscan2_Germline_vcf = filename
                if ".varscan2.somatic.vt.vcf.gz" in filename:
                    if self.TS_varscan2_Somatic_vcf != getime(filename):
                        self.TS_varscan2_Somatic_vcf = getime(filename)
                        self.varscan2_Somatic_vcf = filename
        
        if self.RNA_True == "NA":
            self.RNA_VCF = "NA"
            self.TS_RNA_VCF= "NA"
        else:
            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/alignment/" + self.RNA + "/*.hc.vcf.gz"
            for filename in glob.glob(path):
                if self.TS_RNA_VCF != getime(filename):
                    self.RNA_VCF = filename
                    self.TS_RNA_VCF = getime(filename)
                    self.Gather_RNA_ini()


    def Gather_reports(self):
        if self.DNA_N_True == "NA":
            self.DNA_MultiQC = "NA"
            self.PCGR = "NA"
            self.TS_DNA_MultiQC = "NA"
            self.TS_PCGR = "NA"
        else:
            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna/" + self.Sample + "*multiqc.html"
            for filename in glob.glob(path):
                if self.TS_DNA_MultiQC != getime(filename):
                    self.TS_DNA_MultiQC = getime(filename)
                    self.DNA_MultiQC = filename

            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/ensemble/" + self.Sample + "*/pcgr*/*.flexdb.html"
            for filename in glob.glob(path):
                if self.TS_PCGR != getime(filename):
                    self.TS_PCGR = getime(filename)
                    self.PCGR = filename

        if self.RNA_True == "NA":
            self.RNA_MultiQC = "NA"
            self.TS_RNA_MultiQC = "NA"
        else:
            #####################Not Implemented########################################
            self.RNA_MultiQC = "NA"
            self.TS_RNA_MultiQC = "NA"
            ############################################################################


    def Gather_RNA_other(self):
        if self.RNA_True == "NA":
            self.AnnoFuse = "NA"
            self.GRIDSS= "NA"
            self.RNA_Abundance= "NA"
            self.big_wig_tracks_F = "NA"
            self.big_wig_tracks_R = "NA"
            self.TS_AnnoFuse = "NA"
            self.TS_GRIDSS= "NA"
            self.TS_RNA_Abundance= "NA"
            self.TS_big_wig_tracks_F = "NA"
            self.TS_big_wig_tracks_R = "NA"
        else:
            #####################Not Implemented########################################
            self.GRIDSS= "NA"
            self.TS_GRIDSS= "NA"
            ############################################################################
            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/fusion/" + self.RNA + "/annoFuse/*.putative_driver_fusions.tsv"
            for filename in glob.glob(path):
                if self.TS_AnnoFuse != getime(filename):
                    self.TS_AnnoFuse = getime(filename)
                    self.AnnoFuse = filename

            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/stringtie/" + self.RNA + "/abundance.tab"
            for filename in glob.glob(path):
                if self.TS_RNA_Abundance != getime(filename):
                    self.TS_RNA_Abundance = getime(filename)
                    self.RNA_Abundance = filename

            path = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/tracks/bigWig/" + self.RNA + "*"
            for filename in glob.glob(path):
                if "forward" in filename:
                    if self.big_wig_tracks_F != getime(filename):
                        self.TS_big_wig_tracks_F = getime(filename)
                        self.big_wig_tracks_F = filename
                if "reverse" in filename:
                    if self.big_wig_tracks_R != getime(filename):
                        self.TS_big_wig_tracks_R = getime(filename)
                        self.big_wig_tracks_R = filename

def getime(PATH):
   date = datetime.datetime.fromtimestamp(os.path.getmtime(PATH))
   return date.strftime("%Y/%m/%d")



def main():
    connection = create_connection(r"/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db") 
    ALL_Samples = extract_sample_names(connection)
    for smple in ALL_Samples:
        Sample = Progress(connection,smple)
        update_timestamp_details(Sample) 
        update_fileloc_details(Sample) 
        Sample.Gather_Run_Proc_BAM()
        Sample.Gather_BAM_loc()
        Sample.Gather_Final_BAMs()
        Sample.Gather_VCFs()
        Sample.Gather_reports()
        Sample.Gather_RNA_other()
        update_timestamp_details(Sample)
        update_fileloc_details(Sample)
        Sample.Update_status()
    connection.commit()
    connection.close()


if __name__ == '__main__':
    main()
    #Update db with the objects

