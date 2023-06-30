#!/usr/bin/env bash

variants_bed="$1";
min_pct="$2";
lineage_set="$3";

strains_file="$$_strains.txt"

./filter_variants.awk "$lineage_set" < "$variants_bed" "$strains_file" > "$$_filtered_variants.bed"

sequence_count=$(wc -l "$strains_file" | sed -E "s/^\s*//" | cut -f 1 -d" "); #get just the line count and no leading spaces or file name

# this basically does int(sequence_count*(min_pct/100)) to use min_pct as a percent
scaled_min_count=$(bc <<CALC
scale=10; v=$sequence_count; p=$min_pct/100; s=v*p; scale=0; s=s/1; s;
CALC
)

#echo "$(date +'%b %d %H:%M:%S') - sequence count is $sequence_count, cutoff is $scaled_min_count" >&2

if ((scaled_min_count < 10)); then
  scaled_min_count=10; # exclude variants that occur <10 times in all cases
fi

cut -f 1-3,7,8 "$$_filtered_variants.bed" | uniq -c \
| awk -v min_count="$scaled_min_count" 'BEGIN{ OFS="\t" }; $1 >= min_count { print $2, $3, $4, $5 "/" $6, $1 }';

rm "$$_strains.txt" "$$_filtered_variants.bed";