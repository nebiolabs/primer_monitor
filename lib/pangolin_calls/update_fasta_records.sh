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
\copy pangolin_calls (taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,scorpio_notes,version,pangolin_version,scorpio_version,constellation_version,is_designated,qc_status,qc_notes,note)
  from '$pangolin_csv' WITH (format csv, header on);
UPDATE fasta_records SET fasta_records.$id_field=pangolin_calls.id
  FROM pangolin_calls WHERE pangolin_calls.taxon=fasta_records.strain AND fasta_records.$id_field IS NULL;
CMDS

touch "$1.$(basename "$PWD").complete_pangolin";