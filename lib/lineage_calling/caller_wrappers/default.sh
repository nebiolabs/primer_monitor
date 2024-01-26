#!/usr/bin/env bash

# default no-op lineage caller, ignores input and classifies sequence as "Unknown"

# usage: default.sh <input sequence> [all other arguments ignored]


while getopts ':@:t:T:h' option_arg; do
  case "$option_arg" in
    "@")
      # unused, here to keep a consistent interface for all of these
      ;;
    "t")
      # unused, here to keep a consistent interface for all of these
      ;;
    "T")
      # unused, here to keep a consistent interface for all of these
      ;;
    "h")
      # help
      cat << HELPMSG
Usage: default.sh [-h] <input file>
HELPMSG
      exit 0;
      ;;
    "*")
      echo "invalid option -$option_arg" >&2;
      exit 1;
      ;;
  esac
done

shift $((OPTIND - 1));

if (($# < 1)); then
    echo "usage: default.sh <input file>" >&2;
    exit 1;
fi

input_file="$1"
workdir="$(pwd)";

while read -r line; do
  # name (but not desc) from FASTA header as taxon, lineage of "Unknown", metadata of "default no-op lineage caller"
  echo "$(echo "$line" | cut -f 1 -d " " | sed -E 's/^>//'),Unknown,default no-op lineage caller" >> "$input_file.lineage_calls.csv"
done < <(grep ">" "$workdir/$input_file")
