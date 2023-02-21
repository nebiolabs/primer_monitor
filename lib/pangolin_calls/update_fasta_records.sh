#!/usr/bin/env bash

if (($# < 1)); then
    echo "usage: ./update_fasta_records.sh <pangolin CSV path> [ID field]" >&2;
    exit 1;
fi

if (($# < 2)); then
    id_field="pangolin_call_id";
else
    id_field=$2;
fi

source ../../.env;



psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 <<CMDS;
\copy pangolin_calls from '$1' WITH (format csv, header on);
UPDATE fasta_records SET fasta_records.$id_field=pangolin_calls.id, WHERE pangolin_calls.taxon=fasta_records.strain AND fasta_records.pangolin_call_id IS NULL;
CMDS

touch "$1.$(basename "$PWD").complete_pangolin";