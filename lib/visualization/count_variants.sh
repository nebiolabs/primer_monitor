#!/usr/bin/env bash

variants_bed="$1";
min_count="$2";
threads="$3";

variant_count=$(wc -l < "$variants_bed" | sed -E "s/^\s*//"); #get just the line count and no leading spaces
base_variant_count=100000000;

scaled_min_count=$(((min_count*variant_count)/base_variant_count))

if ((scaled_min_count < 5)); then
  scaled_min_count=5; # exclude variants that occur <5 times in all cases
fi

sort --parallel="$threads" -k2n -k3n < "$variants_bed" | cut -f 1-3,7,8 | uniq -c | awk -v min_count="$scaled_min_count" 'BEGIN{ OFS="\t" }; $1 >= min_count { print $2, $3, $4, $5 "/" $6, $1 }';