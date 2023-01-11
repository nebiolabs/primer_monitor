#!/usr/bin/env bash

mkfifo $$_variants_new; # make named pipe

if (($# == 4)); then # if N variants file provided
    python3 filter_by_date.py "$5" "$3" > $$_variants_new & # filter it by date
elif (($# == 3)); then # otherwise, query for N variants
    ./extract_ns.sh "$3" > $$_variants_new &
else # if anything else, error and print usage
    echo "usage: ns_by_pos.sh <primer BED> <primer set name> <date cutoff> [N variants file]" >&2;
    rm $$_variants_new; # clean up named pipe before exiting
    exit 1;
fi

# make more named pipes
mkfifo $$_ns_cov;

printf "NC_045512.2\t1\t29903\n" > $$_whole_genome.bed; # a single "everything" BED region

bedtools coverage -a $$_whole_genome.bed -b $$_variants_new -d > $$_ns_cov &
python cov_to_bedgraph.py $$_ns_cov;
rm $$_ns_cov $$_variants_new $$_whole_genome.bed; # clean up named pipes and temp files