
import sys
import re

# Usage: ./filter_duplicates.py <old data> <new data> > <output>

old_set = set([])
old_records = open(sys.argv[1], "r")
for line in old_records:
    regex = re.search('\"covv_virus_name\": \"(.*?)\",', line)
    strain = regex.group(1)
    old_set.add(strain)
old_records.close()

new_records = open(sys.argv[2], "r")
for line in new_records:
    regex = re.search('\"covv_virus_name\": \"(.*?)\",', line)
    strain = regex.group(1)
    if strain not in old_set:
        print(line.strip())
old_records.close()
