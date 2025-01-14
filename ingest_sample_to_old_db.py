#!/usr/bin/env python3

import re
import argparse
import csv

from  DB_OPS import create_connection, Update_Samples_Table

def main():
    parser = argparse.ArgumentParser(prog='generate_genpipes_inputs.py', description="Ingest listed samples to old DB.")
    # parser.add_argument(
    #     '--type',
    #     required=True,
    #     help="Type of analysis: either RNA or DNA",
    #     choices=['DNA', 'RNA']
    #     )
    parser.add_argument(
        '-i',
        '--input',
        required=True,
        help="Input file is a csv formatted 'patient,sample_dn,sample_dt,sample_rt'",
        )
    args = parser.parse_args()

    # sequencing_type = args.type

    beluga_db = "/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db"

    #Connect to the db
    connection = create_connection(beluga_db)

    with open(args.input, mode='r') as file:
        # Create a CSV reader object
        csv_reader = csv.DictReader(file)

        # Iterate over each row in the CSV file
        for row in csv_reader:
            patient = row['patient']
            sample_dn = row['sample_dn'] if row['sample_dn'] else 'NA'
            sample_dt = row['sample_dt'] if row['sample_dt'] else 'NA'
            sample_rt = row['sample_rt'] if row['sample_rt'] else 'NA'

            # Determine which sample field to use for parsing
            sample = sample_dn if sample_dn != 'NA' else (sample_dt if sample_dt != 'NA' else sample_rt)

            # Parse cohort and institution
            result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX|HM|CQ)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample)
            cohort = result.group(2)
            institution = result.group(3)

            Update_Samples_Table(
                connection,
                patient,
                patient,
                institution,
                cohort,
                sample_dn,
                sample_dn,
                sample_dt,
                sample_dt,
                sample_rt,
                sample_rt
                )

    connection.commit()
    connection.close()

if __name__ == '__main__':
    main()
