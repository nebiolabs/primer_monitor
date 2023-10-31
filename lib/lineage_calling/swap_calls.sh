#!/usr/bin/env bash

# ensure this is relative to the script's location
source "$(dirname "$0")/../../.env";

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 <<CMDS;
DELETE FROM pangolin_calls WHERE id IN (SELECT pangolin_call_id FROM fasta_records
  WHERE fasta_records.pangolin_call_id IS NOT NULL AND fasta_records.pending_pangolin_call_id IS NOT NULL);
UPDATE fasta_records SET pangolin_call_id=pending_pangolin_call_id, pending_pangolin_call_id=NULL
  WHERE fasta_records.pending_pangolin_call_id IS NOT NULL;
CMDS

touch "$(basename "$PWD").complete_swap";