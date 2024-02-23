#!/usr/bin/env bash

# gets lineages from the DB in CSV format

# you need to export DB_HOST, DB_NAME, and DB_USER_RO before running this

set -e

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -v "slug=$1" <<< "SELECT lineages.name, \
COUNT(fasta_records.id) AS num_seen, min(fasta_records.date_submitted) AS first_seen, \
max(fasta_records.date_submitted) AS last_seen, \
to_timestamp(percentile_cont(0.5) WITHIN GROUP (ORDER BY cast(extract(epoch FROM date_submitted) AS integer)))::date \
AS median_date FROM fasta_records \
INNER JOIN organism_taxa ON organism_taxa.id=fasta_records.organism_taxon_id \
INNER JOIN lineage_calls ON lineage_calls.id=fasta_records.lineage_call_id \
INNER JOIN lineages ON lineages.id=lineage_calls.lineage_id \
INNER JOIN organisms ON organisms.id=organism_taxa.organism_id \
WHERE organisms.slug=:'slug' GROUP BY lineages.id;" --csv -t


