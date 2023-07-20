#!/usr/bin/env bash

# gets BED file of variants from DB, no lineage filtering

set -e

if (($# < 1)); then
    echo "usage: extract_all_variants.sh <date cutoff>" >&2;
    exit 1;
fi
# you need to export DB_HOST, DB_NAME, and DB_USER_RO before running this

# WARNING: do not allow arbitrary user data for date cutoff and lineage set, SQL injection risk

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT organisms.reference_accession, variant_sites.ref_start, variant_sites.ref_end, \
fasta_records.genbank_accession, lineages.name, fasta_records.date_collected, variant_sites.variant_type, \
variant_sites.variant FROM variant_sites \
INNER JOIN fasta_records ON fasta_records.id=variant_sites.fasta_record_id \
INNER JOIN pangolin_calls ON fasta_records.pangolin_call_id=pangolin_calls.id \
INNER JOIN lineages ON pangolin_calls.lineage_id=lineages.id \
INNER JOIN organisms ON organisms.id=lineages.organism_id \
WHERE (date_collected >= '$1
' OR (date_collected IS NULL AND date_submitted >= '$1')) \
ORDER BY variant_sites.ref_start, variant_sites.ref_end, variant_sites.variant_type, variant_sites.variant;" --csv -t \
| tr "," "\t" | grep -v -E "\t[0-9]+N"
