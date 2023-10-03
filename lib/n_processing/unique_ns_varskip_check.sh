#!/usr/bin/env bash

mkfifo $$_variants_new_all_1; # make named pipe
mkfifo $$_variants_new_all_2; # make named pipe

if (($# == 7)); then # if N variants file provided
    tee $$_variants_new_all_1 > $$_variants_new_all_2 < "$7" & # filter it by date
elif (($# == 6)); then # otherwise, query for N variants
    ./extract_ns.sh "$3" | tee $$_variants_new_all_1 > $$_variants_new_all_2 &
else # if anything else, error and print usage
    echo "usage: unique_ns_varskip_check.sh <primer BED> <primer set name> <date cutoff> <N fraction cutoff> <varskip output filename> <non-varskip output filename> [N variants file]" >&2;
    rm $$_variants_new; # clean up named pipe before exiting
    exit 1;
fi

# make more named pipes
mkfifo $$_variants_new_varskip;
mkfifo $$_variants_new_varskip_wc;
mkfifo $$_variants_new_novarskip;

filename="$$_variants_new_varskip_wc"

grep "CDCBI-CRSP" < $$_variants_new_all_1 | tee $$_variants_new_varskip > $filename &


grep -v "CDCBI-CRSP" < $$_variants_new_all_2 | shuf | head -n "$(wc -l $filename | sed -E 's/^[[:space:]]*//g' | cut -f 1 -d ' ')" > $$_variants_new_novarskip &

./process_unique_region_ns.sh "$1" "$2" $$_variants_new_varskip "$4" > "$5";
./process_unique_region_ns.sh "$1" "$2" $$_variants_new_novarskip "$4" > "$6";

rm $$_variants_new_all_1 $$_variants_new_all_2 $$_variants_new_varskip $$_variants_new_varskip_wc $$_variants_new_novarskip; # clean up the named pipes