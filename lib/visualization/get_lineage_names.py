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
        if alias_split[:min(len(search_split), len(alias_split))] == search_split:
            for child in aliases[alias]:
                relevant_aliases += get_relevant_aliases(child, aliases)
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

for relevant in relevant_aliases:
        print(relevant)


