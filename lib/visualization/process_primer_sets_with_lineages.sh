#!/usr/bin/env bash

cutoff_date="$1";
output_path="$2";
min_count="$3";
min_per_primer="$4";
threads="$5";
primer_sets_list_path="$6";
buffer_size="$7";

if [[ $cutoff_date == '-' ]]; then # default of 180 days
  cutoff_date="$(date -d "180 days ago" +"%Y-%m-%d")"
fi

variants_bed="$$_variants.bed"
variants_counts_bed="$$_variants_with_counts.bed";

./extract_all_variants.sh "$cutoff_date" "$threads" "$buffer_size" > "$variants_bed";
shift 7;

for lineage_set_path in "$@"; do
  if [[ $lineage_set_path != "all" ]]; then
    lineages=$(cat "$lineage_set_path");
  else
    lineages="all";
  fi
  lineage_set_name=$(basename "$lineage_set_path" | sed -E "s/\.txt$//")
  echo "processing lineage set $lineage_set_path"
  ./count_variants.sh "$variants_bed" "$min_count" "$threads" "$lineage_set_path" > "$variants_counts_bed";
  xargs ./process_primer_sets.sh "$variants_counts_bed" "$output_path" "$min_per_primer" "$threads" "$lineage_set_name" < "$primer_sets_list_path";
done
