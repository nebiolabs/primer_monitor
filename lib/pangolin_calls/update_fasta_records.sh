#!/usr/bin/env bash

if (($# < 1)); then
    echo "usage: ./update_fasta_records.sh <pangolin CSV path>" >&2;
    exit 1;
fi

psql -d primer_monitor_dev -f update_fasta_records.sql -v pangolin_csv_path=\'"$1"\';
touch "$1.$(basename "$PWD").complete_pangolin";