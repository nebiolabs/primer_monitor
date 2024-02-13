#!/usr/bin/env bash

# nextclade wrapper script

threads=1
temp_dir="$TMPDIR"
dataset_path=""

while getopts ':@:d:t:h' option_arg; do
  case "$option_arg" in
    "@")
      threads="$OPTARG"
      ;;
    "t")
      # temp dir
      temp_dir="$OPTARG"
      ;;
    "d")
      # dataset path
      dataset_path="$OPTARG"
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
workdir="$(pwd)";

if [ -z "$temp_dir" ]; then
  export TMPDIR="$temp_dir"
fi

nextclade run --input-dataset "$BACKEND_INSTALL_PATH/datasets/$dataset_path" \
--output-csv="$input_file.lineage_calls.csv.tmp" \
--jobs "$threads" \
"$workdir/$input_file" # run nextclade

# swap semicolons and commas, remove initial index column
tr ",;" ";," < "$input_file.lineage_calls.csv.tmp" | cut -f 2- -d "," > "$input_file.lineage_calls.csv"
rm "$input_file.lineage_calls.csv.tmp"