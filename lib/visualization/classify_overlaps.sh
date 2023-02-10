#!/usr/bin/env bash

variants_counts_bed="$1";
primer_bed="$2";
min_per_primer="$3";
threads="$4";
intersects_bed="$5";
no_intersects_bed="$6";

bedtools intersect -wa -a "$primer_bed" -b "$variants_counts_bed" | sort --parallel="$threads" | uniq -c \
| awk -v min_per_primer="$min_per_primer" 'BEGIN{ OFS="\t" }; $1 >= min_per_primer { print $2, $3, $4, $5, $1, $7 }' > "$intersects_bed";
bedtools subtract -a "$primer_bed" -b "$intersects_bed" > "$no_intersects_bed";