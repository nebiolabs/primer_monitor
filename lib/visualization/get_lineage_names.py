# Finds all descendant lineages of the given designation
# Usage: curl <pangolin alias_key.json URL> | python get_lineage_names.py <lineage>

# this is also imported as a module by get_lineages_to_show.py

import json
import sys
import re

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

    with open(lineage_filename) as lineages_file:
        has_counts = False
        total_counts = 0
        for line_s in lineages_file:
            if "," in line_s:
                # if the file has a second column, assuming that is counts
                has_counts = True
                line = line_s.split(",")
            else:
                # make this a 1-element list so the rest of the code can stay the same
                line = [line_s.strip()]
            search_str = line[0]+"." # adding trailing "." for same reason as above
            for alias in relevant_aliases_for_search:
                if search_str.startswith(alias):
                    output+=(line[0]+"\n")
                    if has_counts:
                        total_counts += int(line[1])
                    break # break out of the inner loop and go to the next line
    return [output, total_counts]


if __name__ == "__main__":
    result = process_aliases(sys.stdin.read(), sys.argv[2], sys.argv[1])
    print(result[0])
    print("total,"+str(result[1]))
