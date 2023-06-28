# Finds all descendant lineages of the given designation
# Usage: curl <pangolin alias_key.json URL> | python get_lineage_names.py <lineage>

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


aliases_data = json.load(sys.stdin)

search_alias = re.sub('\.?\*$', '', sys.argv[1]) # remove trailing ".*" or "*" from name

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

with open(sys.argv[2]) as lineages_file:
    for line in lineages_file:
        line = line.strip()
        search_str = line+"." # adding trailing "." for same reason as above
        for alias in relevant_aliases_for_search:
            if search_str.startswith(alias):
                print(line)
                break # break out of the inner loop and go to the next line



