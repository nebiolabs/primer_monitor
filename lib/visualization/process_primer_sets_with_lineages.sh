#!/usr/bin/env bash

# entry point

cutoff_date="$1";
output_path="$2";
min_pct="$3";
score_cutoff="$4";
primer_sets_list_path="$5";
organism_path="$6";
threads="$7"

echo "\$DB_NAME is $DB_NAME" >&2

echo "$(date +'%b %d %H:%M:%S') - script started"

mkdir -p "$output_path"

if [[ $cutoff_date == '-' ]]; then # default of 180 days
  cutoff_date="$(date -d "180 days ago" +"%Y-%m-%d")"
fi

variants_bed="$$_variants.bed"
variants_counts_bed="$$_variants_with_counts.bed";
echo "$(date +'%b %d %H:%M:%S') - starting DB fetch"
echo "\"$(dirname "$0")/extract_all_variants.sh\" \"$cutoff_date\" > \"$variants_bed\";"
"$(dirname "$0")/extract_all_variants.sh" "$cutoff_date" > "$variants_bed";
echo "$(date +'%b %d %H:%M:%S') - DB fetch done"
shift 7;

for lineage_set_path in "$@"; do
  lineage_set_name=$(basename "$lineage_set_path" | sed -E "s/\.txt$//")
  echo "$(date +'%b %d %H:%M:%S') - processing lineage set $lineage_set_path"
  echo "\"$(dirname "$0")/count_variants.sh\" \"$variants_bed\" \"$min_pct\" \"$organism_path/lineage_sets/$lineage_set_path\" \"$output_path\" > \"$variants_counts_bed\";" >&2
  "$(dirname "$0")/count_variants.sh" "$variants_bed" "$min_pct" "$organism_path/lineage_sets/$lineage_set_path" "$output_path" > "$variants_counts_bed";
  xargs "$(dirname "$0")/process_primer_sets.sh" "$variants_counts_bed" "$output_path" "$score_cutoff" "$threads" "$lineage_set_name" "$organism_path" < "$primer_sets_list_path";
done
echo "$(date +'%b %d %H:%M:%S') - script done"

rm "$$_variants.bed";