#!/usr/bin/env bash

# ensure this is relative to the script's location
source "$(dirname "$0")/../../.env";

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 <<CMDS;
UPDATE fasta_records SET fasta_records.pangolin_call_id=fasta_records.pending_pangolin_call_id, fasta_records.pending_pangolin_call_id=NULL WHERE fasta_records.pending_pangolin_call_id IS NOT NULL AND fasta_records.pangolin_call_id IS NULL;
CMDS

touch "$(basename "$PWD").complete_swap";