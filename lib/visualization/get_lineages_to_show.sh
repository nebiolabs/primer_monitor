#!/usr/bin/env bash

#Gets "interesting" lineages to display on the igv.js page

#"interesting" = "has at least <cutoff> seqs in it"

#Usage: get_lineages_to_show.py <organism slug> <root lineages> <lineages CSV> <seq counts CSV> <lineage sets path> [overrides path]

# organism slug: the organism slug (to determine which script to run)

#root lineages: a comma-separated string of root lineage names from which to start the tree traversal (e.g. A,B)

#lineages CSV: the CSV of lineage data from get_lineage_data.sh

#seq counts CSV: a file of day-by-day sequence counts for each lineage

#lineage sets path: the path where lineage set definition files (a file listing every Pango lineage recorded in the lineage set) are output

#overrides path: a path to a file listing a set of Pango lineage names to always consider "interesting"


organism_slug="$1"

case "$organism_slug" in
  "sars-cov-2")
    curl -Ssf https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json \
      | python "$(dirname "$0")/per_organism_scripts/get_lineages_to_show_sars_cov_2.py" A,B "$@"
    ;;
  "rsv")
    cat <(curl -Ssf https://raw.githubusercontent.com/rsv-lineages/lineage-designation-A/main/.auto-generated/clades.md) \
        <(curl -Ssf https://raw.githubusercontent.com/rsv-lineages/lineage-designation-B/main/.auto-generated/clades.md) \
        | python "$(dirname "$0")/per_organism_scripts/get_lineages_to_show_rsv.py" A,B "$@"
    ;;
  *)
    python "$(dirname "$0")/per_organism_scripts/get_lineages_to_show_default.py" "" "$@"
    ;;
esac


