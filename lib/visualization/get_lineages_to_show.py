# Gets "interesting" lineages to display on the igvjs page
# "interesting" = "has at least <cutoff> seqs in it"

import get_lineage_names
import sys

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
def get_lineages(alias_str, lineages_file, starting_lineages, cutoff, lineages = set(), root_lineages = None):
    # if this is the initial call
    if root_lineages is None:
        # set the root_lineages list so neither base clade is dropped
        root_lineages = set(starting_lineages)
    for lineage in starting_lineages:
        # get total seqs in this lineage group, and all lineage names
        result, total = get_lineage_names.process_aliases(alias_str, lineages_file, lineage)
        # if this lineage is "interesting" (has more than cutoff # of seqs or is a root lineage)
        if (total>cutoff or lineage in root_lineages) and lineage not in lineages:
            # add it and repeat recursively on that lineage
            lineages.add(lineage)
            lineages = get_lineages(alias_str, lineages_file, get_sublineages(result, lineage), cutoff, lineages, root_lineages)
    return lineages

starting_lineages = sys.argv[1].split(",")
alias_str = sys.stdin.read()
lineages_file = sys.argv[2]
cutoff = int(sys.argv[3])

interesting_lineages = get_lineages(alias_str, lineages_file, starting_lineages, cutoff)

for lineage in interesting_lineages:
    print(lineage)

