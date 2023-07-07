#!/usr/bin/env bash

variants_counts_bed="$1";
output_path="$2";
score_cutoff="$3";
threads="$4";
name="$5";
organism_path="$6";

shift 6;

for primer_bed in "$@"; do
	outputname=$(basename "$primer_bed" | sed -E "s/\.bed$//I");
	mkdir -p "$output_path/primer_sets/${outputname}";
	"$(dirname "$0")/primers_affected.sh" "$variants_counts_bed" "$organism_path/primer_sets_raw/$primer_bed" "$output_path/primer_sets/${outputname}/$name.bed" "$score_cutoff" "$threads" ;
done

rm "$variants_counts_bed";