#!/usr/bin/env bash

# pangolin wrapper script

#  <taxon record ID> is unused and is present to keep the interface of this and nextclade compatible

threads=1
temp_dir="$TMPDIR"

while getopts ':@:d:t:h' option_arg; do
  case "$option_arg" in
    "@")
      threads="$OPTARG"
      ;;
    "t")
      # temp dir
      temp_dir="$OPTARG"
      ;;
    "h")
      # help
      cat << HELPMSG
Usage: pangolin.sh [-@ THREADS] [-t TEMP_DIR] [-h] <input file>
HELPMSG
      exit 0;
      ;;
    "d")
      # unused, here to keep a consistent interface for all of these
      ;;
    "*")
      echo "invalid option -$option_arg" >&2;
      exit 1;
      ;;
  esac
done

shift $((OPTIND - 1));

if (($# < 1)); then
    echo "usage: pangolin.sh [-@ THREADS] [-t TEMP_DIR] [-h] <input file>" >&2;
    exit 1;
fi

input_file="$1"
threads="$2"
temp_dir="$3"
workdir="$(pwd)";

pangolin "$workdir/$input_file" -t "$threads" -o "$workdir" \
--outfile "$input_file.lineage_calls.csv" \
${temp_dir:+"--tempdir $temp_dir"}; # run pangolin
