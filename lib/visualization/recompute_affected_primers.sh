#!/usr/bin/env bash

# Actually recompute what primers are affected (the outer script mostly does prep stuff and uploads the results from this one)

# Usage: recompute_affected_primers.sh <cutoff date> <output path> <frequency cutoff> <score cutoff> <primer sets list> <organism dir path> <cpus> [lineage sets...]

# if cutoff date is "-", default to 180 days

set -e

cutoff_date="$1";
output_path="$2";
min_pct="$3";
score_cutoff="$4";
primer_sets_list_path="$5";
organism_slug="$6";
threads="$7"

# defines echo_log function
source "$(dirname "$0")/../echo_log.sh"

echo_log "overlap computation started"

mkdir -p "$output_path"

if [[ $cutoff_date == '-' ]]; then # default to lookback date from DB
  lookback=$(psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -v "organism_slug=$organism_slug" \
  <<< "SELECT variant_bed_lookback_days FROM organisms WHERE slug=:'organism_slug';" --csv -t)
  cutoff_date="$(date -d "$lookback days ago" +"%Y-%m-%d")"
fi

variants_bed="$$_variants.bed"
variants_counts_bed="$$_variants_with_counts.bed";
echo_log "starting DB fetch"
"$(dirname "$0")/extract_all_variants.sh" "$cutoff_date" "$organism_slug" > "$variants_bed";
echo_log "DB fetch done"
shift 7;

for lineage_set_path in "$@"; do
  lineage_set_name=$(basename "$lineage_set_path" | sed -E "s/\.txt$//")
  echo_log "processing lineage set $lineage_set_path ($lineage_set_name)"
  "$(dirname "$0")/count_variants.sh" "$variants_bed" "$min_pct" "./$organism_slug/lineage_sets/$lineage_set_path" "$output_path" > "${lineage_set_name}_$variants_counts_bed";
  xargs "$(dirname "$0")/process_primer_sets.sh" "${lineage_set_name}_$variants_counts_bed" "$output_path" "$score_cutoff" "$threads" "$lineage_set_name" "./$organism_slug" < "$primer_sets_list_path";
done
echo_log "overlap computation done"

rm "$$_variants.bed";