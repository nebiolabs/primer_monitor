"""
Gets "interesting" lineages to display on the igv.js page

default version, returns everything (which should just be "Unknown" or similar)

Usage: get_lineages_to_show.py <root lineages> <lineages CSV> <seq counts CSV> <lineage sets path> [overrides path]

root lineages: (ignored here)

lineages CSV: the CSV of lineage data from get_lineage_data.sh. Returns all lineages in this file.

seq counts CSV: (ignored here)

lineage sets path: the path where lineage set definition files (a file listing every lineage recorded in the lineage set) are output

overrides path: a path to a file listing a set of Pango lineage names to always consider "interesting"
"""

import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('root_lineages')
parser.add_argument('lineages_csv')
parser.add_argument('seq_counts_csv')
parser.add_argument('output_path')
parser.add_argument('organism_slug')
parser.add_argument('overrides_path', nargs='?')
parsed_args = parser.parse_args()

# ignores the other args
lineage_list_file = parsed_args.lineages_csv

output_path = parsed_args.output_path

overrides = []

if parsed_args.overrides_path is not None:
    overrides_file = parsed_args.overrides_path
    with open(overrides_file) as f:
        for override_line in f:
            overrides.append(override_line.strip())

interesting_lineages = set()

with open(lineage_list_file) as f:
    for line_s in f:
        lineage_rec = line_s.strip().split(",")
        interesting_lineages.add(lineage_rec[0])

for override_lineage in overrides:
    interesting_lineages.add(override_lineage)

lineage_groups = {}

print("{")
for lineage_group in interesting_lineages:
    recorded_lineages = lineage_group
    display_name = lineage_group
    with open(output_path + "/" + lineage_group + ".txt", "w") as f:
        f.write(recorded_lineages)
    print('"' + lineage_group + '": "' + display_name + '",')
print('"all": "All"')
print("}")
