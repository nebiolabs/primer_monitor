#!/usr/bin/env bash

cutoff_date=$1;
output_path=$2;
min_count=$3;
min_per_primer=$4;
threads=$5;

variants_bed="$$_variants.bed"
variants_counts_bed="$$_variants_with_counts.bed";

shift 5;

./extract_all_variants.sh "$cutoff_date" > "$variants_bed";

./count_variants.sh "$variants_bed" "$variants_counts_bed";

rm "$variants_bed";

for primer_bed in "$@"; do
	outputname=$(echo "$primer_bed" | sed -E "s/\.bed$//I");
	mkdir -p "$output_path/${outputname}";
	./primers_affected.sh "$variants_counts_bed" "$primer_bed" "$output_path/${outputname}/color.bed" "$min_per_primer" "$threads";
done

rm "$variants_counts_bed";
