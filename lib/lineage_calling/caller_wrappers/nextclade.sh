#!/usr/bin/env bash

# nextclade wrapper script

threads=1
temp_dir="$TMPDIR"
dataset=""

while getopts ':@:t:d:h' option_arg; do
  case "$option_arg" in
    "@")
      threads="$OPTARG"
      ;;
    "t")
      # a list of primer sets to process
      temp_dir="$OPTARG"
      ;;
    "d")
      # the name of the nextclade dataset
      dataset="$OPTARG"
      ;;
    "h")
      # help
      cat << HELPMSG
Usage: nextclade.sh [-@ THREADS] [-t TEMP_DIR] [-d DATASET] [-h] <input file>
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
    echo "usage: nextclade.sh [-@ THREADS] [-t TEMP_DIR] [-d DATASET] [-h] <input file>" >&2;
    exit 1;
fi

input_file="$1"
threads="$2"
temp_dir="$3"
workdir="$(pwd)";

pangolin "$workdir/$input_file" -t "$threads" -o "$workdir" \
--outfile "$input_file.lineage_calls.csv" \
${temp_dir:+"--tempdir $temp_dir"}; # run pangolin