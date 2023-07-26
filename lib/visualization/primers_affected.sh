#!/usr/bin/env bash

set -e

variants_counts_bed="$1";
primer_bed="$2";
output_file="$3";
score_cutoff="$4";
threads="$5";

"$(dirname "$0")/classify_overlaps.sh" "$variants_counts_bed" "$primer_bed" "$score_cutoff" "$threads" "$$_affected.bed" "$$_unaffected.bed";

awk 'BEGIN{ OFS="\t" }; { print $1, $2, $3, $4, $5, $6, $2, $2, "0,0,200" }' < "$$_unaffected.bed" > "$$_unaffected_color.bed";
awk 'BEGIN{ OFS="\t" }; { print $1, $2, $3, $4, $5, $6, $2, $2, "230,100,0" }' < "$$_affected.bed" > "$$_affected_color.bed";

sort -k1 -k2n "$$_unaffected_color.bed" "$$_affected_color.bed" > "$output_file";

rm "$$_unaffected.bed" "$$_affected.bed" "$$_unaffected_color.bed" "$$_affected_color.bed";
