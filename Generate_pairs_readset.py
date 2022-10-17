#!/usr/bin/env python3
import glob
import sys
import os
import re
import datetime
import shutil
from  DB_OPS import create_connection, Update_Samples_Table

def main():
    try:
        Mol_type = sys.argv[1]
    except IndexError:
        print ("You must specify DNA or RNA")
        sys.exit(1)
    Mol_type = sys.argv[1]
    print (Mol_type)
    if Mol_type != 'DNA' and Mol_type != 'RNA':
        raise NameError('You must specify DNA or RNA')
    
    #file locations.
    Output_RR = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN/raw_reads/"
    Work_dir = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/"
    Input_RR = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/raw_reads/"
    SQL_DB = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db"
    #Output_RR = "/home/dipop/MOH/TEST/MAIN/raw_reads/"
    #Work_dir = "/home/dipop/MOH/TEST/MAIN/"
    #Input_RR = "/home/dipop/MOH/TEST/raw_reads/"
    #SQL_DB = "/home/dipop/MOH/TEST/TEST.db"
    date = datetime.date.today()
    TIME = date.strftime("%Y-%m-%d")   
    #Connect to the db
    connection = create_connection(SQL_DB)

    #Populate the lists
    RNA_Samples = []
    Tumour_Samples = []
    Normal_Samples = []
    Bad_Names = []
    All_files = os.listdir(Input_RR)
    for name in All_files:
        #Test the name. Throw warning if file is bad. Toss it in Bad names.
        Tester = re.match("^MoHQ-(JG|CM|GC|MU|MR|XX)-\w+-\w+-\w+-\w+(D|R)(T|N)", name)
        if Tester == None:
            print (f"{name} is in inproper format")
            Bad_Names.append(name)
        #Test if it is RNA
        elif (name[-2] == "R"):
            RNA_Samples.append(name)
        #Test if it is DNA and a Tumour
        elif (name[-2] == "D" and name[-1] == "T" ):
            Tumour_Samples.append(name)
        #Test if it is DNA and Normal
        elif (name[-2] == "D" and name[-1] == "N" ):
            Normal_Samples.append(name)

    #Make file for bad names
    if Bad_Names:
        with open(Work_dir + "bad_names.txt", "w+") as f:
            for line in Bad_Names:
                f.write(f"{line}\n")
        print (f"We found poorly named files in {Input_RR}. A list of files that need to be corrected/reprocessed can be found at {Work_dir}bad_names.txt")
    #RNA 
    if Mol_type == 'RNA':
        #Check to see if any files have already present in the final directory
        duplicates = []
        for name in RNA_Samples:
            if os.path.exists(Output_RR + name):
                duplicates.append(name)
        if duplicates:
            with open(Work_dir + "duplicates.txt", "w+") as f:
                for line in duplicates:
                    f.write(f"{line}\n")
            print (f"We found directories in both your input {Input_RR} and the output {Output_RR} A list of files that need to be corrected/reprocessed can be found at {Work_dir}duplicates.txt. ")
            sys.exit("Program will not produce readsets until this is solved. Exiting") 

        #readset Generation
        if RNA_Samples:
            readset = []
            readset.append('Sample	Readset	LibraryType	RunType	Run	Lane	Adapter1	Adapter2	QualityOffset	BED	FASTQ1	FASTQ2	BAM\n')
            for name in RNA_Samples:
                mylines = []                             
                with open (f'{Input_RR}{name}/{name}_readset.tsv', 'rt') as myfile: 
                    for myline in myfile:                
                        mylines.append(myline)             
                readset = readset + mylines[1:]
            with open(f"{Work_dir}{TIME}_RNA_readset.tsv", "w+") as f:
                for line in readset:
                    f.write(f"{line}")
                print (f"Generated {Work_dir}{TIME}_RNA_readset.tsv")
            
            #move files
            for name in RNA_Samples:
                shutil.move(f"{Input_RR}{name}",f"{Output_RR}{name}")

            #Add To database
            for name in RNA_Samples:
                result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", name)
                Sample = result.group(1)
                Cohort = result.group(2)
                Institution = result.group(3)
                RNA_Samp = name
                Update_Samples_Table(connection,Sample,Sample,Institution,Cohort,"NA","NA","NA","NA",RNA_Samp,RNA_Samp)

        else:
            print("No RNA Samples to Move")

    #DNA
    elif Mol_type == 'DNA':
        #Check to see if any files have already present in the final directory
        DNA_Samples = Tumour_Samples + Normal_Samples
        duplicates = []
        for name in DNA_Samples:
            if os.path.exists(Output_RR + name):
                duplicates.append(name)
        if duplicates:
            with open(Work_dir + "duplicates.txt", "w+") as f:
                for line in duplicates:
                    f.write(f"{line}\n")
            print (f"We found directories in both your input {Input_RR} and the output {Output_RR} A list of files that need to be corrected/reprocessed can be found at {Work_dir}duplicates.txt. ")
            sys.exit("Program will not produce readsets or pairs until this is solved. Exiting") 
        
        #Find matching pairs
        pair_dict = {}
        #print("Normal len: " + str(len(Normal_Samples)))
        #print("Tumour len: " + str(len(Tumour_Samples)))
        for Normal in Normal_Samples:
            #print("Normal: " + Normal)
            result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", Normal)
            Sample = result.group(1)
            for Tumour in Tumour_Samples:
                result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", Tumour)
                patient = result.group(1)
                #print("Tumour: " + Tumour)
                if Sample == patient:
                    print("Pair found for sample " + Sample)
                    pair_dict[Normal] = Tumour
        DNA_Samples = [] 
        if pair_dict:
            readset = []
            readset.append('Sample	Readset	LibraryType	RunType	Run	Lane	Adapter1	Adapter2	QualityOffset	BED	FASTQ1	FASTQ2	BAM\n')
            pairs = []
            for normal in pair_dict:
                DNA_Samples.append(normal)
                tumor = pair_dict[normal]
                DNA_Samples.append(tumor)
                #readset
                mylines = []                             
                with open (f'{Input_RR}{normal}/{normal}_readset.tsv', 'rt') as myfile: 
                    for myline in myfile:                
                        mylines.append(myline)             
                readset = readset + mylines[1:]
                mylines = []                             
                with open (f'{Input_RR}{tumor}/{tumor}_readset.tsv', 'rt') as myfile: 
                    for myline in myfile:                
                        mylines.append(myline)             
                readset = readset + mylines[1:]
               #Add to db 
                result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", normal)
                Sample = result.group(1)
                Cohort = result.group(2)
                Institution = result.group(3)
                Update_Samples_Table(connection,Sample,Sample,Institution,Cohort,normal,normal,tumor,tumor,"NA","NA")

            #pairs
                pairs.append(f"{Sample},{normal},{tumor}")

            #move file
                shutil.move(f"{Input_RR}{normal}",f"{Output_RR}{normal}")
                shutil.move(f"{Input_RR}{tumor}",f"{Output_RR}{tumor}")

            #write the files
            with open(f"{Work_dir}{TIME}_TP_readset.tsv", "w+") as f:
                for line in readset:
                    f.write(f"{line}")
                print (f"Generated {Work_dir}{TIME}_TP_readset.tsv")
            with open(f"{Work_dir}{TIME}_TP_pairs.csv", "w+") as f:
                for line in pairs:
                    f.write(f"{line}\n")
                print (f"Generated {Work_dir}{TIME}_TP_pairs.csv")

        else:
            print("No DNA pairs to Move")


    connection.commit()
    connection.close()

if __name__ == '__main__':
    main()
