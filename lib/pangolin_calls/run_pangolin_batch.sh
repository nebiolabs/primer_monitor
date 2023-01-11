#!/usr/bin/env bash

if (($# < 6)); then
    echo "usage: ./run_pangolin.sh <input path> <output file name> <done file path> <output dir path> <pangolin path> <threads>" >&2;
    exit 1;
fi

PATH=$PATH:$5/bin "$5/bin/pangolin" "$1" -t "$6" -o "$4" --outfile "$2"; # run pangolin
touch "$3"; # create file to indicate this is done
