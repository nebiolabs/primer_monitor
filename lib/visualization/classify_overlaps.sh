#!/usr/bin/env bash

set -e

# Classifies primers as affected/unaffected for a single primer set

# Usage: classify_overlaps.sh <variant counts bed> <primers bed> <score cutoff> <threads> <affected output filename> <unaffected output filename>

variants_counts_bed="$1";
primer_bed="$2";
score_cutoff="$3";
threads="$4";
filter_pass_bed="$5";
filter_fail_bed="$6";

bedtools intersect -wa -wb -a "$primer_bed" -b "$variants_counts_bed" | sort --parallel="$threads" > "$$_raw_intersects.bed"

python3 "$(dirname "$0")/classify_overlaps.py" "$$_raw_intersects.bed" "$filter_pass_bed" "$$_filter_fail_intersects.bed" "$score_cutoff";

bedtools subtract -a "$primer_bed" -b "$$_raw_intersects.bed" > "$$_no_intersects.bed";

cat "$$_no_intersects.bed" "$$_filter_fail_intersects.bed" | sort -k2n > "$filter_fail_bed";

rm "$$_raw_intersects.bed" "$$_filter_fail_intersects.bed" "$$_no_intersects.bed";