#!/usr/bin/env bash

# you need to export DB_HOST, DB_NAME, and DB_USER_RO before running this

set -e

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT lineages.name, \
COUNT(fasta_records.id) AS num_seen, min(fasta_records.date_submitted) AS first_seen, \
max(fasta_records.date_submitted) AS last_seen, \
to_timestamp(percentile_cont(0.5) WITHIN GROUP (ORDER BY cast(extract(epoch FROM date_submitted) AS integer)))::date \
AS median_date FROM fasta_records \
INNER JOIN pangolin_calls ON pangolin_calls.id=fasta_records.pangolin_call_id \
INNER JOIN lineages ON lineages.id=pangolin_calls.lineage_id GROUP BY lineages.id;" --csv -t


