#!/usr/bin/env bash

# counts occurrences of variants

set -e

variants_bed="$1";
min_pct="$2";
lineage_set="$3";
output_path="$4";

strains_file="$$_strains.txt"

#just in case there are no strains, make an empty file
touch "$strains_file"
lineage_name=$(basename "$lineage_set" | sed -E "s/\.txt$//");



# handle "all" lineage set properly
if [ "$lineage_name" == "all" ]; then
  lineage_set="all"
fi

"$(dirname "$0")/filter_variants.awk" "$lineage_set" < "$variants_bed" "$(pwd)/$strains_file" > "$$_filtered_variants.bed"


mkdir -p "$output_path/lineage_variants"
LINECOUNT=$(wc -l < "$$_filtered_variants.bed");
awk -F"\t" '{ OFS="\t"; print $1, $2, $3, $7 "->" $8 }' "$$_filtered_variants.bed" | sort -k2n -k3n -k4 | uniq -c \
| awk -v linecount="$LINECOUNT" '$1 > linecount*0.001 { OFS="\t"; print $2, $3, $4, $5 }' > "$output_path/lineage_variants/$lineage_name.bed"


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