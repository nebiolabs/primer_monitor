#!/usr/bin/env bash

if (($# < 2)); then
    echo "usage: ./run_pangolin.sh <input file name> <threads> [temp dir]" >&2;
    exit 1;
fi

workdir=$(pwd);
pangolin "$workdir/$1" -t "$2" -o "$workdir" --outfile "$1.pangolin_calls.csv" ${3:+"--tempdir"} "$3"; # run pangolin
touch "$workdir/$1.$(basename "$PWD").pangolin_calls.done"; # create file to indicate this is done