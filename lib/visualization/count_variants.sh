#!/usr/bin/env bash

variants_bed="$1";
min_count="$2";
threads="$3";
lineage_set="$4";

./filter_variants.awk "$lineage_set" < "$variants_bed" > "$$_filtered_variants.bed"

variant_count=$(wc -l < "$$_filtered_variants.bed" | sed -E "s/^\s*//"); #get just the line count and no leading spaces
base_variant_count=100000000;

scaled_min_count=$(((min_count*variant_count)/base_variant_count))

#echo "$(date +'%b %d %H:%M:%S') - variant count is $variant_count, cutoff is $scaled_min_count" >&2

if ((scaled_min_count < 10)); then
  scaled_min_count=10; # exclude variants that occur <10 times in all cases
fi

cut -f 1-3,7,8 "$$_filtered_variants.bed" | uniq -c \
| awk -v min_count="$scaled_min_count" 'BEGIN{ OFS="\t" }; $1 >= min_count { print $2, $3, $4, $5 "/" $6, $1 }';