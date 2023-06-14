#!/usr/bin/env bash

if (($# < 1)); then
    pending="";
else
    pending=$1;
fi

source ../../.env;

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 <<CMDS;
UPDATE fasta_records SET fasta_records.pangolin_call_id=fasta_records.pending_pangolin_call_id, fasta_records.pending_pangolin_call_id=NULL, WHERE fasta_records.pending_pangolin_call_id IS NOT NULL AND fasta_records.pangolin_call_id IS $pending NULL;
CMDS

touch "$(basename "$PWD").complete_swap";