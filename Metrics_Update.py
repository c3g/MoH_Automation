#!/usr/bin/env python3
import glob
import sys
import re
from  DB_OPS import update_metrics_db,create_connection,extract_sample_field,extract_sample_names
#TO DO::
#FIX!

def main():
    connection = create_connection(r"/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db") 
    
    Samples = extract_sample_names(connection)
    #TEST CASE
    #Samples = ['MoHQ-JG-9-10','MoHQ-MU-12-122','MoHQ-MU-12-1224']

    All_Samples = []
    PAIRED_SAMPLES = dict()
    for sample in Samples:
        All_Samples.append(extract_sample_field(connection,sample,"DNA_N"))
        PAIRED_SAMPLES[extract_sample_field(connection,sample,"DNA_N")] = sample
        All_Samples.append(extract_sample_field(connection,sample,"DNA_T"))
        PAIRED_SAMPLES[extract_sample_field(connection,sample,"DNA_T")] = sample
        All_Samples.append(extract_sample_field(connection,sample,"RNA"))
        PAIRED_SAMPLES[extract_sample_field(connection,sample,"RNA")] = sample
    All_Samples = [x for x in All_Samples if x != ""]
    extract_data(All_Samples,connection,PAIRED_SAMPLES)
    connection.commit()
    connection.close()

def extract_data(SAMP,connection,PAIRED_SAMPLES):
    for Sample in SAMP:
        PATIENT = PAIRED_SAMPLES[Sample]
        TYPE=''
        if Sample.endswith('N'):
            TYPE = 'N'
        elif Sample.endswith('T'):
            TYPE = 'T'
        Y_Flags=[]
        R_Flags=[]
        WGS_Bases_Over_Q30 = extract_bs_over_q30(Sample,PATIENT)
        if WGS_Bases_Over_Q30 == None:
            WGS_Bases_Over_Q30 = ''
        elif (float(WGS_Bases_Over_Q30)<75):
            R_Flags.append('WGS Bases Over Q30')
        elif (float(WGS_Bases_Over_Q30)<80):
            Y_Flags.append('WGS Bases Over Q30')

        WGS_Min_Aligned_Reads_Delivered = extract_min_aln_rds(Sample,PATIENT)
        if WGS_Min_Aligned_Reads_Delivered == None or WGS_Min_Aligned_Reads_Delivered == '' :
            WGS_Min_Aligned_Reads_Delivered = ''            
        elif (float(WGS_Min_Aligned_Reads_Delivered)<260000000 and TYPE == 'N'):
            R_Flags.append('WGS Min Aligned Reads Delivered')
        elif (float(WGS_Min_Aligned_Reads_Delivered)<660000000 and TYPE == 'N'):
            Y_Flags.append('WGS Min Aligned Reads Delivered')
        elif (float(WGS_Min_Aligned_Reads_Delivered)<530000000 and TYPE == 'T'):
            R_Flags.append('WGS Min Aligned Reads Delivered')
        elif (float(WGS_Min_Aligned_Reads_Delivered)<1330000000 and TYPE == 'T'):
            Y_Flags.append('WGS Min Aligned Reads Delivered')
        

        #WGS_Duplication_Rate
        WGS_duplicates = extract_sambama_dups(Sample,PATIENT)
        if  WGS_duplicates == None or WGS_duplicates == '' :
            WGS_duplicates = ''            
        elif (float(WGS_duplicates)>50):
            R_Flags.append('WGS_Duplication_Rate')
        elif (float(WGS_duplicates)>20):
            Y_Flags.append('WGS_Duplication_Rate')

        WGS_Raw_Coverage = extract_raw_coverage(Sample)
        if  WGS_Raw_Coverage == None or WGS_Raw_Coverage == '' :
            WGS_Raw_Coverage = ''            
        elif (float(WGS_Raw_Coverage)<30 and TYPE == 'N'):
            R_Flags.append('WGS Raw Coverage')
        elif (float(WGS_Raw_Coverage)<80 and TYPE == 'T'):
            R_Flags.append('WGS Raw Coverage')


        WGS_Dedup_Coverage = extract_ded_coverage(Sample,PATIENT)
        if WGS_Dedup_Coverage  == None or WGS_Dedup_Coverage == '' :
            WGS_Dedup_Coverage = ''            

        Median_Insert_Size = extract_insert_size(Sample,PATIENT)
        if  Median_Insert_Size == None or Median_Insert_Size == '' :
            Median_Insert_Size = ''            
        elif (float(Median_Insert_Size)<300):
            Y_Flags.append('Median_Insert_Size')
        elif (float(Median_Insert_Size)<150):
            R_Flags.append('Median_Insert_Size')

        #WGS_Contamination
        WGS_Contamination = extract_contamination(Sample,PATIENT)
        if WGS_Contamination  == None or WGS_Contamination == '' :
            WGS_Contamination = ''            
        elif (float(WGS_Contamination)>5):
            R_Flags.append('WGS_Contamination')

        #Concordance
        Concordance = extract_concordance(Sample,PATIENT)
        if  Concordance == None or Concordance == '' :
            Concordance = ''            
        elif (float(Concordance)<99):
            R_Flags.append('Concordance')

        #Tumor_Purity
        Purity = extract_purity(Sample,PATIENT)
        if  Purity == None or Purity == '' :
            Purity = ''            
        elif (float(Purity)<30):
            R_Flags.append('Purity')

        #WTS_Clusters
        WTS_Clusters = extract_WTS_Clusters(Sample)
        if  WTS_Clusters == None or WTS_Clusters == '' :
           WTS_Clusters  = ''            
        elif (float(WTS_Clusters)<80000000):
            R_Flags.append('WTS_Clusters')
        elif (float(WTS_Clusters)<100000000):
            Y_Flags.append('WTS_Clusters')

        #WTS_Exonic_Rate
        WTS_Exonic_Rate = extract_WTS_exonic(Sample)
        if WTS_Exonic_Rate  == None or WTS_Exonic_Rate == '' :
            WTS_Exonic_Rate = ''            
        elif (float(WTS_Exonic_Rate)<0.6):
            R_Flags.append('WTS_Exonic_Rate')
        elif (float(WTS_Exonic_Rate)<0.8):
            Y_Flags.append('WTS_Exonic_Rate')

        #WTS_Unique_Reads
        WTS_Unique_Reads = extract_WTS_unique(Sample)
        if  WTS_Unique_Reads == None or WTS_Unique_Reads == '' :
            WTS_Unique_Reads = ''            

        #WTS_rRA_contamination
        rRNA_count = extract_WTS_rRNA(Sample)
        if  rRNA_count == None or rRNA_count == '' or WTS_Unique_Reads == '':
            WTS_rRNA_contamination = ''            
        else:
            WTS_rRNA_contamination = int(rRNA_count)/int(WTS_Unique_Reads)
            if (float(WTS_rRNA_contamination)>0.35):
                R_Flags.append('WTS_rRNA_contamination')
            elif (float(WTS_rRNA_contamination)>0.1):
                Y_Flags.append('WTS_rRNA_contamination')

        mean_ins_size = extract_mean_map_ins_size(Sample)
        med_ins_size = extract_med_map_ins_size(Sample)

        #Warning flags
        Yellow_Flags=';'.join(Y_Flags)
        Red_Flags=';'.join(R_Flags)
        
        update_metrics_db(connection,Sample,WGS_Bases_Over_Q30,WGS_Min_Aligned_Reads_Delivered,WGS_Raw_Coverage,WGS_Dedup_Coverage,Median_Insert_Size,WGS_duplicates,WGS_Contamination,WTS_Clusters,WTS_Unique_Reads,WTS_Exonic_Rate,WTS_rRNA_contamination,Concordance,Purity,Yellow_Flags,Red_Flags,mean_ins_size,med_ins_size)

        print (Sample)
        #print (WGS_Bases_Over_Q30)
        #print (WGS_Min_Aligned_Reads_Delivered)
        #print (WGS_Raw_Coverage)
        #print (WGS_Dedup_Coverage)
        #print (WGS_duplicates)
        print (Median_Insert_Size)
        print (mean_ins_size)
        print (med_ins_size)
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



def extract_WTS_rRNA(ID):
    if 'DT' in ID or 'DN' in ID or 'D-T' in ID or 'D-N' in ID:
        return ''
    else:
        path = '/home/dipop/MOH/MAIN/metrics/rna/' + ID + '/rnaseqc/*/*rRNA_counts.txt'
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                lines=f.readlines()
                line=lines[0]
                fields = line.split("\t")
                return fields[0] 

def extract_WTS_unique(ID):
    if 'DT' in ID or 'DN' in ID or 'D-T' in ID or 'D-N' in ID:
        return ''
    else:
        path = '/home/dipop/MOH/MAIN/metrics/rna/' + ID + '/rnaseqc/*/*.metrics.tmp.txt'
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                lines=f.readlines()
                line=lines[3]
                fields = line.split("\t")
                Output = fields[0]
                return Output
    
def extract_WTS_exonic(ID):
    if 'DT' in ID or 'DN' in ID or 'D-T' in ID or 'D-N' in ID:
        return ''
    else:
        path = '/home/dipop/MOH/MAIN/metrics/rna/' + ID + '/rnaseqc/*/*.metrics.tmp.txt'
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                lines=f.readlines()
                line=lines[7]
                fields = line.split("\t")
                Output = float(fields[1])*100
                return (f"%.2f" % round(Output, 2))

def extract_WTS_Clusters(ID):
    if 'DT' in ID or 'DN' in ID or 'D-T' in ID or 'D-N' in ID:
        return ''
    else:
        path = '/home/dipop/MOH/MAIN/metrics/run_metrics/*'
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                for line in f:
                    if ID in line:
                        data = line.split(",")
                        return data[12]

def extract_raw_coverage(ID):
    if 'R' in ID:
        return ''
    else:
        path = '/home/dipop/MOH/MAIN/metrics/run_metrics/*'
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                for line in f:
                    if ID in line:
                        data = line.split(",")
                        return data[41]
def extract_med_map_ins_size(ID):
    path = '/home/dipop/MOH/MAIN/metrics/run_metrics/*'
    for filename in glob.glob(path):
        with open(filename, 'r') as f:
            for line in f:
                if ID in line:
                    data = line.split(",")
                    return data[38]
def extract_mean_map_ins_size(ID):
    path = '/home/dipop/MOH/MAIN/metrics/run_metrics/*'
    for filename in glob.glob(path):
        with open(filename, 'r') as f:
            for line in f:
                if ID in line:
                    data = line.split(",")
                    return data[39]

def extract_purity(ID,PATIENT):
    if 'R' in ID or ID.endswith('N'):
        return ''
    else:
        path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/pairedVariants/' + PATIENT + '/purple/*.purity.tsv'    
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                lines=f.readlines()
                line=lines[1]
                fields = line.split("\t")
                Output = float(fields[0])*100
                return Output 

#WGS N % Contamination,WGS T % Contamination,
def extract_contamination(ID,PATIENT):
    if 'R' in ID:
        return ''
    else:
        path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/' + PATIENT + '*.contamination.tsv'    
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                for line in f:
                    words = line.split(" ")
                    output = words[-1]
                    output = output[:-2]
                    if words[0].startswith('N') and ID.endswith('N'):
                        return output
                    elif words[0].startswith('T') and ID.endswith('T'):
                        return output


def extract_concordance(ID,PATIENT):
    if 'R' in ID:
        return ''
    else:
        CON = None
        path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/' + PATIENT + '*.concordance.tsv'    
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                CON = f.readline()
        if CON:
            CON = CON.replace('Concordance:','')
            CON = CON.replace(' ','')
            CON = CON[:-1]
            CON = CON[:-1]
            if float(CON)<=1:
                CON=float(CON)*100
            return CON
        else:
            return None


def extract_sambama_dups(ID,PATIENT):
    if 'R' in ID:
        return ''
    else:
        DUPS = None;
        path = '/home/dipop/MOH/MAIN/metrics/run_metrics/*'
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                for line in f:
                    if ID in line:
                        data = line.split(",")
                        DUPS = data[15]
        if DUPS == None or DUPS == '':
            return ''
        else:
            return float(DUPS) 

def extract_insert_size(ID,PATIENT):
    if 'R' in ID:
        OUTPUT = 0;
        path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna/' +  ID + '/*.insert_size_metrics'
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                lines=f.readlines()
                line=lines[7]
                metrics = line.split("\t")
                OUTPUT = metrics[0]
        return OUTPUT 
    else:
        return parse_multiqc2(ID,'median_insert_size',2,PATIENT)


def extract_ded_coverage(ID,PATIENT):
    if 'R' in ID:
        return ''
    else:
        path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna/' +  ID + '/qualimap/' + ID + '/genome_results.txt'
        OUTPUT = ''
        for filename in glob.glob(path):
            with open(filename, 'r') as f:
                lines=f.readlines()
                line=lines[71]
                #line is      mean coverageData = 304.9902X
                metrics = line.split(" ")
                OUTPUT = metrics[-1].replace('X', '')
        if OUTPUT == '' or OUTPUT == '':
            return ""
        else:
            return float(OUTPUT)

def extract_min_aln_rds(ID,PATIENT):
    if 'R' in ID:
        return ''
    else:
        return parse_multiqc(ID,'QualiMap_mqc-generalstats-qualimap-mapped_reads',0,PATIENT)

def parse_multiqc(ID,field,Round,PATIENT):
    path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna/' +  PATIENT + '*.multiqc_data/multiqc_general_stats.txt'
    for filename in glob.glob(path):
        with open(filename, 'r') as f:
            Output =''
            header = f.readline().split("\t")
            index = header.index(field)
            A = f.readline().split("\t")
            B = f.readline().split("\t")
            if A[0].endswith(ID[-1]):
                Output = A[index]
            elif B[0].endswith(ID[-1]):
                Output = B[index]
            if Output != '':
                Output=float(Output)
                Output=f"%.{Round}f" % round(Output, Round)
            return Output

def parse_multiqc2(ID,field,Round,PATIENT):
    path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna/' +  PATIENT + '*.multiqc_data/multiqc_qualimap_bamqc_genome_results.txt'
    for filename in glob.glob(path):
        with open(filename, 'r') as f:
            Output =''
            header = f.readline().split("\t")
            index = header.index(field)
            A = f.readline().split("\t")
            B = f.readline().split("\t")
            if A[0].endswith(ID[-1]):
                Output = A[index]
            elif B[0].endswith(ID[-1]):
                Output = B[index]
            if Output != '':
                Output=float(Output)
                Output=f"%.{Round}f" % round(Output, Round)
            return Output

def extract_bs_over_q30(ID,PATIENT):
    path = ''
    if 'R' in ID:
        path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/rna/' + ID + '/picard_metrics/*quality_distribution_metrics'
    else:
        path = '/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/metrics/dna/' + ID + '*/picard_metrics/*quality_distribution_metrics'
    Tester = re.compile('(\d+)\W+(\d+)')
    for filename in glob.glob(path):
        with open(filename, 'r') as f:
            Abv_30=0
            Blw_30=0
            for line in f:
                if line[:1].isdigit():
                    Test = Tester.match(line)
                    Qual = Test.group(1)
                    count = Test.group(2)
                    if (int(Qual) < 30):
                        Blw_30 += int(count)
                    else :
                        Abv_30 += int(count)
            percent_abv=(Abv_30/(Abv_30+Blw_30))*100
            percent_abv = "%.2f" % round(percent_abv, 2)
            return percent_abv




if __name__ == '__main__':
    main()
    #Update db with the objects

