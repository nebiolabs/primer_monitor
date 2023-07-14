#!/usr/bin/env bash

set -e

variants_counts_bed="$1";
primer_bed="$2";
score_cutoff="$3";
threads="$4";
intersects_bed="$5";
no_intersects_bed="$6";

bedtools intersect -wa -wb -a "$primer_bed" -b "$variants_counts_bed" | sort --parallel="$threads" > "$$_raw_intersects.bed"

python3 "$(dirname "$0")/classify_overlaps.py" "$$_raw_intersects.bed" "$intersects_bed" "$score_cutoff";

bedtools subtract -a "$primer_bed" -b "$intersects_bed" > "$no_intersects_bed";

rm "$$_raw_intersects.bed"