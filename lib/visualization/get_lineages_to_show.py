# Gets "interesting" lineages to display on the igvjs page
# "interesting" = "has at least <cutoff> seqs in it"

import get_lineage_names
import sys
import json
import datetime

def passes_filter(lineage, total, cutoff, lineage_data):
    today = datetime.datetime.now().date()
    if lineage_data[lineage]['last_seen'] is None:
        # no recorded lineages, fail
        return False

    timediff = (today-lineage_data[lineage]['median_date']).days
    return total > cutoff*(timediff/90)

# gets all sublineages for a given lineage listed in the lineages file
def get_sublineages(result, base):
    base_depth = len(base.split("."))
    prefixes = set()
    # loop over lineage names
    for rec in result.split("\n"):
        line = rec.strip()
        # if this is a child of the starting lineage
        if line.startswith(base+"."):
            # add 1 level deeper than base
            prefix = ".".join(line.split(".")[:base_depth+1])
            prefixes.add(prefix)
        # if this is something else (e.g. an alias)
        elif not base.startswith(line):
            # add the prefix
            prefix = line.split(".")[0]
            prefixes.add(prefix)
        else:
            # don't add anything
            pass
    return prefixes

# recursively traverses the lineages tree to find "interesting" lineages
def get_lineages(alias_str, lineages_file, starting_lineages, cutoff, lineage_data, lineages = set(), root_lineages = None):
    # if this is the initial call
    if root_lineages is None:
        # set the root_lineages list so neither base clade is dropped
        root_lineages = set(starting_lineages)
    for lineage in starting_lineages:
        # get total seqs in this lineage group, and all lineage names
        result, total, min_date, max_date, median_date = get_lineage_names.process_aliases(alias_str, lineages_file, lineage)
        # if this lineage is "interesting" (has more than cutoff # of seqs or is a root lineage)
        if lineage not in lineage_data:
            lineage_data[lineage] = {'num_seen':total,
                                     'first_seen':min_date,
                                     'last_seen':max_date,
                                     'median_date':median_date}
        if (passes_filter(lineage, total, cutoff, lineage_data) or lineage in root_lineages) and lineage not in lineages:
            # add it and repeat recursively on that lineage
            lineages.add(lineage)
            lineages = get_lineages(alias_str, lineages_file, get_sublineages(result, lineage), cutoff, lineage_data, lineages, root_lineages)
    return lineages

starting_lineages = sys.argv[1].split(",")
alias_str = sys.stdin.read()
lineages_file = sys.argv[2]
cutoff = int(sys.argv[3])
output_path = sys.argv[4]
lineage_data = {}
with open(sys.argv[2]) as f:
    for line_s in f:
        line = line_s.strip().split(",")
        lineage_data[line[0]] = {'num_seen':line[1],
                                 'first_seen':datetime.datetime.strptime(line[2], '%Y-%m-%d').date(),
                                 'last_seen':datetime.datetime.strptime(line[3], '%Y-%m-%d').date(),
                                 'median_date':datetime.datetime.strptime(line[4], '%Y-%m-%d').date()}

interesting_lineages = get_lineages(alias_str, lineages_file, starting_lineages, cutoff, lineage_data)

lineage_groups = {}

print("{")
for lineage in interesting_lineages:
    lineages = get_lineage_names.process_aliases(alias_str, lineages_file, lineage)[0]
    with open(output_path+"/"+lineage+".txt", "w") as f:
        f.write(lineages)
    print('"'+lineage+'*": "'+lineage+'",')
print("}")

