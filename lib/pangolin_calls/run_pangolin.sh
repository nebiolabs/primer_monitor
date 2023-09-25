#!/usr/bin/env bash

if (($# < 2)); then
    echo "usage: ./run_pangolin.sh <input file name> <threads> [temp dir]" >&2;
    exit 1;
fi

input_file="$1"
threads="$2"
temp_dir="$3"
workdir="$(pwd)";

pangolin "$workdir/$input_file" -t "$threads" -o "$workdir" \
--outfile "$input_file.pangolin_calls.csv" \
${temp_dir:+"--tempdir $temp_dir"}; # run pangolin

touch "$workdir/$input_file.$(basename "$PWD").pangolin_calls.done"; # create file to indicate this is done