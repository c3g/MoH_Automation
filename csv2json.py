#!/usr/bin/env python3

import os
import csv
import json
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert CSV from pgweb to JSON for modification')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    parser.add_argument('-o', '--output', required=False, help='Output JSON file')

    return parser.parse_args()

def convert_csv_to_json(input_file, output_file):
    # Initialize the dictionary with the specified format
    data = {
        "modification": []
    }

    # Read the CSV file and process the headers
    with open(input_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)

        # Extract the string before '_' in the headers and add to the dictionary
        for header in headers:
            table_name = header.split('_')[0]
            if table_name not in [mod["table"] for mod in data["modification"]]:
                data["modification"].append({"table": table_name, "id": []})

        # Fill the id field with the content of the CSV file
        for row in reader:
            for i, value in enumerate(row):
                if value:  # Ensure no id is empty
                    table_name = headers[i].split('_')[0]
                    for mod in data["modification"]:
                        if mod["table"] == table_name:
                            mod["id"].append(value)

    # Remove duplicate ids from lists
    for mod in data["modification"]:
        mod["id"] = list(set(mod["id"]))

    # Convert the dictionary to a JSON string
    json_data = json.dumps(data, indent=4)

    # Write the JSON string to the output file
    with open(output_file, 'w') as jsonfile:
        jsonfile.write(json_data)

def main():
    args = parse_arguments()

    output = args.output or f"{os.path.basename(args.input)}".replace('.csv', '.json')

    convert_csv_to_json(args.input, output)

if __name__ == "__main__":
    main()
