#!/usr/bin/env bash

# Process overlaps for primer sets for a single lineage set

# Usage: process_primer_sets.sh <variant counts bed> <output path> <score cutoff> <threads> <lineage set name> <organism dir path> [primer sets...]

set -e

variants_counts_bed="$1";
output_path="$2";
score_cutoff="$3";
threads="$4";
name="$5";
organism_path="$6";

shift 6;

for primer_bed in "$@"; do
  outputname=$(basename "$primer_bed" | sed -E "s/\.bed$//I");
  mkdir -p "$output_path/primer_sets_status/${outputname}";
  "$(dirname "$0")/primers_affected.sh" "$variants_counts_bed" "$organism_path/primer_sets_bed/$primer_bed" "$output_path/primer_sets_status/${outputname}/$name.bed" "$score_cutoff" "$threads";
done

rm "$variants_counts_bed";