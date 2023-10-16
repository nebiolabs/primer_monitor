#!/usr/bin/env bash

# ensure this is relative to the script's location
source "$(dirname "$0")/../../.env";

# deletes any pangolin call not connected to a fasta record

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 <<CMDS;
DELETE FROM pangolin_calls WHERE NOT EXISTS
  (SELECT 1 FROM fasta_records WHERE fasta_records.pangolin_call_id=pangolin_calls.id
  OR fasta_records.pending_pangolin_call_id=pangolin_calls.id);
CMDS