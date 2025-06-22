#!/usr/bin/env python
import argparse
import json
import os

# Parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Merge results from multiple temp result files")
    parser.add_argument('-d', '--directory', type=str, default="Null", help="directory of the simulation")
    parser.add_argument('-n', '--num_files', type=int, default=1, help="number of temp files to merge")
    args = parser.parse_args()
    if args.directory[-1] != '/':
        args.directory += '/'
    return args

Args = parse_args()

print(f"Merging ({Args.num_files}) temp results...")
# Initialize an empty dictionary to hold the merged results
merged_results = {}

# Loop through each file and merge its contents
for i in range(Args.num_files):
    with open(f"{Args.directory}temp_titration_result_{i+1}.json", 'r') as file:
        temp_data = json.load(file)
        merged_results.update(temp_data)

# Save the merged results to titration_result.json
with open(f"{Args.directory}titration_result.json", 'w') as file:
    json.dump(merged_results, file)

for i in range(Args.num_files):
    os.remove(f"{Args.directory}temp_titration_result_{i+1}.json")

print("Merged results saved to titration_result.json")