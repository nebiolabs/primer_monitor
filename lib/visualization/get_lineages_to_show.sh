#!/usr/bin/env bash

#Gets "interesting" lineages to display on the igv.js page

#"interesting" = "has at least <cutoff> seqs in it"

#Usage: get_lineages_to_show.py <caller name> <root lineages> <lineages CSV> <seq counts CSV> <lineage sets path> [overrides path]

# caller name: the name of the lineage caller used here

#root lineages: a comma-separated string of root lineage names from which to start the tree traversal (e.g. A,B)

#lineages CSV: the CSV of lineage data from get_lineage_data.sh

#seq counts CSV: a file of day-by-day sequence counts for each lineage

#lineage sets path: the path where lineage set definition files (a file listing every Pango lineage recorded in the lineage set) are output

#overrides path: a path to a file listing a set of Pango lineage names to always consider "interesting"


caller_name="$1"
shift;

case "$caller_name" in
  pangolin)
    curl -Ssf https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json \
      | python "$(dirname "$0")/per_caller_scripts/pangolin/get_lineages_to_show.py" "$@"
    ;;
  *)
    python "$(dirname "$0")/per_caller_scripts/default/get_lineages_to_show.py" "$@"
    ;;
esac


