#!/usr/bin/env bash

cutoff_date="$1";
output_path="$2";
min_count="$3";
min_per_primer="$4";
threads="$5";
primer_sets_list_path="$6";

shift 6;

for lineage_set_path in "$@"; do
  lineage_set_name=$(basename "$lineage_set_path" | sed -E "s/\.csv$//")
  echo "processing lineage set $lineage_set_path"
  xargs ./process_primer_sets.sh "$cutoff_date" "$output_path" "$min_count" "$min_per_primer" "$threads" "$lineage_set_name" "$lineage_set_path" < "$primer_sets_list_path";
done
