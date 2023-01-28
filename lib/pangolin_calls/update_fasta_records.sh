#!/usr/bin/env bash

if (($# < 1)); then
    echo "usage: ./update_fasta_records.sh <pangolin CSV path>" >&2;
    exit 1;
fi

source ../../.env;

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 <<CMDS;
CREATE TEMPORARY TABLE pangolin_tmp (taxon varchar NOT NULL, lineage varchar NOT NULL, conflict varchar, ambiguity_score float4, scorpio_call varchar, scorpio_support numeric, scorpio_conflict numeric, scorpio_notes text, version varchar NOT NULL, pangolin_version varchar NOT NULL, scorpio_version varchar NOT NULL, constellation_version varchar NOT NULL, is_designated text NOT NULL, qc_status varchar, qc_notes varchar, note varchar);
\copy pangolin_tmp from '$1' WITH (format csv, header on);
UPDATE fasta_records SET pangolin_lineage=pangolin_tmp.lineage, pangolin_version=pangolin_tmp.pangolin_version FROM pangolin_tmp WHERE pangolin_tmp.taxon=fasta_records.strain;
DROP TABLE pangolin_tmp;
CMDS

touch "$1.$(basename "$PWD").complete_pangolin";