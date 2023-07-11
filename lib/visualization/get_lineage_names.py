# Finds all descendant lineages of the given designation
# Usage: curl <pangolin alias_key.json URL> | python get_lineage_names.py <lineage>

# this is also imported as a module by get_lineages_to_show.py

import json
import sys
import re
import statistics
import datetime

def get_relevant_aliases(search_string, aliases):
    relevant_aliases = [search_string]
    for alias in aliases.keys():
        search_split = search_string.split(".")
        alias_split = alias.split(".")
        search_len = len(search_split)
        alias_len = len(alias_split)
        min_len = min(search_len, alias_len)
        if alias_split[:min_len] == search_split[:min_len]:
            extra_segment = ""
            # if search string is more specific than alias
            if search_len > alias_len:
                # add the extra specificity onto the alias
                extra_segment = ".".join([""]+search_split[alias_len:search_len])
            for child in aliases[alias]:
                relevant_aliases += get_relevant_aliases(child+extra_segment, aliases)
    return relevant_aliases

def process_aliases(alias_str, lineage_filename, search_string):

    output = ""

    aliases_data = json.loads(alias_str)

    search_alias = re.sub('\.?\*$', '', search_string) # remove trailing ".*" or "*" from name

    reversed_aliases = {}

    for alias in aliases_data:
        if type(aliases_data[alias]) == type([]): # is a list
            for parent in aliases_data[alias]:
                if parent not in reversed_aliases:
                    reversed_aliases[parent] = []
                reversed_aliases[parent].append(alias)
        else:
            parent = aliases_data[alias]
            if parent not in reversed_aliases:
                reversed_aliases[parent] = []
            reversed_aliases[parent].append(alias)


    relevant_aliases = get_relevant_aliases(search_alias, reversed_aliases)

    relevant_aliases_for_search = []

    for relevant in relevant_aliases:
        relevant_aliases_for_search.append(relevant+".") # adding trailing "." so e.g. "A" doesn't match "AY"

    lineage_data = [[], [], [], []]
    with open(lineage_filename) as lineages_file:
        has_counts = False
        has_dates = False
        for line_s in lineages_file:
            if "," in line_s:
                # if the file has a second column
                line = line_s.split(",")
                #print(line)
                #print(len(line))
                if len(line) >= 2:
                    has_counts = True
                if len(line) >= 5:
                    has_dates = True
            else:
                # make this a 1-element list so the rest of the code can stay the same
                line = [line_s.strip()]
            search_str = line[0]+"." # adding trailing "." for same reason as above
            for alias in relevant_aliases_for_search:
                #print(alias, search_str)
                if search_str.startswith(alias):
                    output+=(line[0]+"\n")
                    #print("checking "+str(line))
                    if has_counts:
                        lineage_data[0].append(int(line[1]))
                    if has_dates:
                        lineage_data[1].append(datetime.datetime.strptime(line[2], '%Y-%m-%d').date())
                        lineage_data[2].append(datetime.datetime.strptime(line[3], '%Y-%m-%d').date())
                        lineage_data[3].append(datetime.datetime.strptime(line[4].strip(), '%Y-%m-%d').date())
                    break # break out of the inner loop and go to the next line
    min_date = min(lineage_data[1]) if len(lineage_data[1]) > 0 else None
    max_date = max(lineage_data[2]) if len(lineage_data[2]) > 0 else None
    median_date = datetime.datetime.fromtimestamp(statistics.median([datetime.datetime(year=date_rec.year,month=date_rec.month,day=date_rec.day).timestamp() for date_rec in lineage_data[3]])).date() if len(lineage_data[3]) > 0 else None
    return [output, sum(lineage_data[0]), min_date, max_date, median_date]


if __name__ == "__main__":
    result = process_aliases(sys.stdin.read(), sys.argv[2], sys.argv[1])
    print(result[0])
    print("total,"+str(result[1]))
