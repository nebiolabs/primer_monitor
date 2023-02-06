#!/usr/bin/env bash

variants_counts_bed="$1";
primer_bed="$2";
output_file="$3";
min_per_primer="$4";
threads="$5";

bedtools intersect -wa -a "$primer_bed" -b "$variants_counts_bed" | sort --parallel="$threads" | uniq -c | awk -v min_per_primer="$min_per_primer" 'BEGIN{ OFS="\t" }; $1 >= min_per_primer { print $2, $3, $4, $5, $1, $7 }' > $$_intersects.bed;
bedtools subtract -a "$primer_bed" -b $$_intersects.bed > $$_no_intersects.bed;
awk 'BEGIN{ OFS="\t" }; { print $1, $2, $3, $4, $5, $6, $2, $2, "0,0,200" }' < $$_no_intersects.bed > $$_no_intersects_color.bed;
awk 'BEGIN{ OFS="\t" }; { print $1, $2, $3, $4, $5, $6, $2, $2, "230,100,0" }' < $$_intersects.bed > $$_intersects_color.bed;
sort -k1 -k2n $$_no_intersects_color.bed $$_intersects_color.bed > "$output_file";

rm $$_no_intersects.bed $$_intersects.bed $$_no_intersects_color.bed $$_intersects_color.bed;