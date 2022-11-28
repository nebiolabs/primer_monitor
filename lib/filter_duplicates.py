
import sys
import re

# Usage: ./filter_duplicates.py <old data> <new data> > <output>

def get_strain(line):
    regex = re.search('\"covv_virus_name\": \"(.*?)\",', line)
    return regex.group(1)

record_set = set([])
with open(sys.argv[1], "r") as old_records:
    for line in old_records:
        if line.strip() != "": # skip blank lines
            strain = get_strain(line)
            record_set.add(strain)

with open(sys.argv[2], "r") as new_records:
    for line in new_records:
        line_stripped = line.strip() # so line.strip() isn't called twice
        if line_stripped != "": # skip blank lines
            strain = get_strain(line)
            if strain not in record_set:
                print(line_stripped)
                record_set.add(strain) # do not output two of the same record, old or new
