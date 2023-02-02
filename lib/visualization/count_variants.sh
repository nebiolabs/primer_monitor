#!/usr/bin/env bash

variants_bed=$1;
output_file=$2;
min_count=$3;

sort --parallel="$threads" -k2n -k3n < "$variants_bed" | cut -f 1-3,6,7 | uniq -c | awk -v min_count="$min_count" 'BEGIN{ OFS="\t" }; $1 >= min_count { print $2, $3, $4, $5 "/" $6, $1 }' > "$output_file";