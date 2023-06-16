#!/usr/bin/env bash

if (($# < 1)); then
    echo "usage: ./update_fasta_records.sh <pangolin CSV path> [ID field]" >&2;
    exit 1;
fi

pangolin_csv=$1;

if (($# < 2)); then
    id_field="pangolin_call_id";
else
    id_field=$2;
fi

# ensure this is relative to the script's location
source "$(dirname "$0")/../../.env";



psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 <<CMDS;
\copy pangolin_calls from '$pangolin_csv' WITH (format csv, header on);
UPDATE fasta_records SET fasta_records.$id_field=pangolin_calls.id, WHERE pangolin_calls.taxon=fasta_records.strain AND fasta_records.$id_field IS NULL;
CMDS

touch "$1.$(basename "$PWD").complete_pangolin";