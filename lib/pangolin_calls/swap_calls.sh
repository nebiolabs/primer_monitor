#!/usr/bin/env bash

if (($# < 1)); then
    echo "usage: ./swap_calls.sh [pending_only]" >&2;
    exit 1;
fi

if (($# < 2)); then
    pending="";
else
    pending=$2;
fi

source ../../.env;

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 <<CMDS;
UPDATE fasta_records SET fasta_records.pangolin_call_id=fasta_records.pending_pangolin_call_id, fasta_records.pending_pangolin_call_id=NULL, WHERE fasta_records.pending_pangolin_call_id IS NOT NULL AND fasta_records.pangolin_call_id IS $pending NULL;
CMDS

touch "$1.$(basename "$PWD").complete_swap";