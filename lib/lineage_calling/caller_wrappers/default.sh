#!/usr/bin/env bash

# default no-op lineage caller, ignores input and classifies sequence as "Unknown"

# usage: default.sh <input sequence> [all other arguments ignored]

if (($# < 1)); then
    echo "usage: default.sh <input file name>" >&2;
    exit 1;
fi

input_file="$1"
workdir="$(pwd)";

while read -r line; do
  # name (but not desc) from FASTA header as taxon, lineage of "Unknown", metadata of "default no-op lineage caller"
  echo "$(echo "$line" | cut -f 1 -d " "),Unknown,default no-op lineage caller" >> "$input_file.lineage_calls.csv"
done < <(grep ">" "$workdir/$input_file")
