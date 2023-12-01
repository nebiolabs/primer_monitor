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

# ignores the other args
lineage_list_file = sys.argv[2]

output_path = sys.argv[4]

overrides = []

if len(sys.argv) >= 6:
    overrides_file = sys.argv[5]
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
