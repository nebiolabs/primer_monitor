#!/usr/bin/env bash

# recompute the primer data for igvjs visualization
set -e # fail on error

primer_monitor_path="$1"
organism_dirname="$2"
pct_cutoff="$3"
score_cutoff="$4"
cpus="$5"

# unset any pre-existing value for JUMP_PROXY so unset == "don't use a jump proxy"
unset JUMP_PROXY

source "$primer_monitor_path/.env"
export DB_HOST
export DB_USER_RO
export DB_NAME
mkdir -p "$organism_dirname/config"
mkdir -p "$organism_dirname/lineage_sets"
mkdir -p "$organism_dirname/lineage_variants"
mkdir -p "$organism_dirname/primer_sets_raw"

"$primer_monitor_path/lib/visualization/get_lineage_data.sh" > lineages.csv

"$primer_monitor_path/lib/visualization/get_primer_sets.sh" "$organism_dirname/primer_sets_raw" > "$organism_dirname/config/tracks.json"

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT COALESCE(date_collected, date_submitted), COUNT(*) \
FROM fasta_records GROUP BY COALESCE(date_collected, date_submitted);" --csv -t > seq_counts.csv

curl https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json \
| python "$primer_monitor_path/lib/visualization/get_lineages_to_show.py" A,B lineages.csv seq_counts.csv "$organism_dirname/lineage_sets"

cat <(printf "{") <(find "$organism_dirname/lineage_sets" -exec basename -a "{}" + \
| sed -E 's/^(.*)\\.txt$/"\\1": "\\1.*",/') <(echo '"all": "All"}') > "$organism_dirname/config/lineage_sets.json"

ls "$organism_dirname/primer_sets_raw" > primer_sets_data.txt

cat <(ls "$organism_dirname/lineage_sets") <(echo "all.txt") \
| xargs "$primer_monitor_path/lib/visualization/process_primer_sets_with_lineages.sh" - "./$organism_dirname" "$pct_cutoff" "$score_cutoff" \
primer_sets_data.txt "./$organism_dirname" "$cpus"

rm primer_sets_data.txt

# remove old files so this doesn't clutter up the directories

# shellcheck disable=SC2029
# I actually want the variables expanded on the client side
ssh ${JUMP_PROXY:+"-J"} "${JUMP_PROXY:-''}"  "$FRONTEND_HOST" "rm -rf $IGVSTATIC_PATH/$organism_dirname/primer_sets; \
rm -f $IGVSTATIC_PATH/$organism_dirname/primer_sets_raw/* $IGVSTATIC_PATH/$organism_dirname/lineage_sets/* \
$IGVSTATIC_PATH/$organism_dirname/lineage_variants/*;"

# copies over the new files
scp ${JUMP_PROXY:+"-o"} ${JUMP_PROXY:+"ProxyJump=$JUMP_PROXY"} -r "./$organism_dirname/*" "$FRONTEND_HOST:$IGVSTATIC_PATH/$organism_dirname/";