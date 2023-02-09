#!/usr/bin/env bash

variants_counts_bed="$1";
output_path="$2";
min_per_primer="$3";
threads="$4";
name="$5";

shift 5;

for primer_bed in "$@"; do
	outputname=$(echo "$primer_bed" | sed -E "s/\.bed$//I");
	mkdir -p "$output_path/${outputname}";
	./primers_affected.sh "$variants_counts_bed" "$primer_bed" "$output_path/${outputname}/$name.bed" "$min_per_primer" "$threads" ;
done

rm "$variants_counts_bed";
