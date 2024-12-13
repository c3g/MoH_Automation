#!/usr/bin/env python3

import sys
import os
import re
import datetime
import shutil
from  DB_OPS import create_connection, Update_Samples_Table

def main():
    try:
        seq_type = sys.argv[1]
    except IndexError:
        print ("You must specify DNA or RNA")
        sys.exit(1)
    seq_type = sys.argv[1]
    print (seq_type)
    if seq_type != 'DNA' and seq_type != 'RNA':
        raise NameError('You must specify DNA or RNA')

    #file locations.
    moh_processing_folder = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING"
    main_raw_reads_folder = os.path.join(moh_processing_folder, "MAIN", "raw_reads")
    transferred_raw_reads_folder = os.path.join(moh_processing_folder, "raw_reads")
    # SQL_DB = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db"
    #main_raw_reads_folder = "/home/dipop/MOH/TEST/MAIN/raw_reads/"
    #moh_processing_folder = "/home/dipop/MOH/TEST/MAIN/"
    #transferred_raw_reads_folder = "/home/dipop/MOH/TEST/raw_reads/"
    #SQL_DB = "/home/dipop/MOH/TEST/TEST.db"
    date = datetime.date.today()
    time = date.strftime("%Y-%m-%d")
    #Connect to the db
    connection = create_connection("/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db")

    #Populate the lists
    rna_samples = []
    tumour_samples = []
    normal_samples = []
    bad_names = []
    all_files = os.listdir(transferred_raw_reads_folder)
    for name in all_files:
        #Test the name. Throw warning if file is bad. Toss it in Bad names.
        tester = re.match("^MoHQ-(JG|CM|GC|MU|MR|XX|CQ)-\w+-\w+-\w+-\w+(D|R)(T|N)", name)
        if tester == None:
            print (f"{name} is in inproper format")
            bad_names.append(name)
        #Test if it is RNA
        elif name[-2] == "R":
            rna_samples.append(name)
        #Test if it is DNA and a tumour
        elif (name[-2] == "D" and name[-1] == "T" ):
            tumour_samples.append(name)
        #Test if it is DNA and normal
        elif (name[-2] == "D" and name[-1] == "N" ):
            normal_samples.append(name)

    #Make file for bad names
    if bad_names:
        with open(os.path.join(moh_processing_folder, "bad_names.txt"), "w+") as file:
            for line in bad_names:
                file.write(f"{line}\n")
        print (f"We found poorly named files in {transferred_raw_reads_folder}. A list of files that need to be corrected/reprocessed can be found at {moh_processing_folder}bad_names.txt")
    #RNA
    if seq_type == 'RNA':
        #Check to see if any files have already present in the final directory
        duplicates = []
        for name in rna_samples:
            if os.path.exists(os.path.join(main_raw_reads_folder, name)):
                duplicates.append(name)
        if duplicates:
            with open(os.path.join(moh_processing_folder, "duplicates.txt"), "w+") as file:
                for line in duplicates:
                    file.write(f"{line}\n")
            print (f"We found directories in both your input {transferred_raw_reads_folder} and the output {main_raw_reads_folder} A list of files that need to be corrected/reprocessed can be found at {moh_processing_folder}duplicates.txt. ")
            sys.exit("Program will not produce readsets until this is solved. Exiting")

        #readset Generation
        if rna_samples:
            readset = []
            readset.append('Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM\n')
            for name in rna_samples:
                mylines = []
                with open (os.path.join(transferred_raw_reads_folder, name, name + '_readset.tsv'), 'rt') as myfile:
                    for myline in myfile:
                        mylines.append(myline)
                readset = readset + mylines[1:]
            rna_readset = os.path.join(moh_processing_folder, time + "_RNA_readset.tsv")
            with open(rna_readset, "w+") as file:
                for line in readset:
                    file.write(f"{line}")
                print (f"Generated {rna_readset}")

            #move files
            for name in rna_samples:
                shutil.move(os.path.join(transferred_raw_reads_folder, name), os.path.join(main_raw_reads_folder, name))

            #Add To database
            for name in rna_samples:
                result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX|CQ)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", name)
                sample = result.group(1)
                cohort = result.group(2)
                institution = result.group(3)
                rna_samp = name
                Update_Samples_Table(connection,sample,sample,institution,cohort,"NA","NA","NA","NA",rna_samp,rna_samp)

        else:
            print("No RNA Samples to Move")

    #DNA
    elif seq_type == 'DNA':
        #Check to see if any files have already present in the final directory
        dna_samples = tumour_samples + normal_samples
        duplicates = []
        for name in dna_samples:
            if os.path.exists(os.path.join(main_raw_reads_folder, name)):
                duplicates.append(name)
        if duplicates:
            duplicates_file = os.path.join(moh_processing_folder, "duplicates.txt")
            with open(duplicates_file, "w+") as f:
                for line in duplicates:
                    f.write(f"{line}\n")
            print (f"We found directories in both your input {transferred_raw_reads_folder} and the output {main_raw_reads_folder} A list of files that need to be corrected/reprocessed can be found at {duplicates_file}.")
            sys.exit("Program will not produce readsets or pairs until this is solved. Exiting")

        #Find matching pairs
        pair_dict = {}
        #print("normal len: " + str(len(normal_samples)))
        #print("tumour len: " + str(len(tumour_samples)))
        for normal in normal_samples:
            #print("normal: " + normal)
            result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX|CQ)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", normal)
            sample = result.group(1)
            for tumour in tumour_samples:
                result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX|CQ)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", tumour)
                patient = result.group(1)
                #print("tumour: " + tumour)
                if sample == patient:
                    print("Pair found for sample " + sample)
                    pair_dict[normal] = tumour
        dna_samples = []
        if pair_dict:
            readset = []
            readset.append('Sample\tReadset\tLibraryType\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM\n')
            pairs = []
            for normal in pair_dict:
                dna_samples.append(normal)
                tumor = pair_dict[normal]
                dna_samples.append(tumor)
                #readset
                mylines = []
                with open (os.path.join(transferred_raw_reads_folder, normal, normal + '_readset.tsv'), 'rt') as myfile:
                    for myline in myfile:
                        mylines.append(myline)
                readset = readset + mylines[1:]
                mylines = []
                with open (os.path.join(transferred_raw_reads_folder, tumor, tumor + '_readset.tsv'), 'rt') as myfile:
                    for myline in myfile:
                        mylines.append(myline)
                readset = readset + mylines[1:]
               #Add to db
                result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX|CQ)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", normal)
                sample = result.group(1)
                cohort = result.group(2)
                institution = result.group(3)
                Update_Samples_Table(connection,sample,sample,institution,cohort,normal,normal,tumor,tumor,"NA","NA")

            #pairs
                pairs.append(f"{sample},{normal},{tumor}")

            #move file
                shutil.move(os.path.join(transferred_raw_reads_folder, normal), os.path.join(main_raw_reads_folder, normal))
                shutil.move(os.path.join(transferred_raw_reads_folder, tumor), os.path.join(main_raw_reads_folder, tumor))

            #write the files
            tp_readset_file = os.path.join(moh_processing_folder, time + "_TP_readset.tsv")
            tp_pair_file = os.path.join(moh_processing_folder, time + "_TP_pairs.csv")
            with open(tp_readset_file, "w+") as file:
                for line in readset:
                    file.write(f"{line}")
                print (f"Generated {tp_readset_file}")
            with open(tp_pair_file, "w+") as file:
                for line in pairs:
                    file.write(f"{line}\n")
                print (f"Generated {tp_pair_file}")

        else:
            print("No DNA pairs to Move")


    connection.commit()
    connection.close()

if __name__ == '__main__':
    main()
