#!/usr/bin/env bash

# gets BED file of variants from DB, no lineage filtering

set -e

if (($# < 1)); then
    echo "usage: extract_all_variants.sh <date cutoff> <organism slug>" >&2;
    exit 1;
fi
# you need to export DB_HOST, DB_NAME, and DB_USER_RO before running this

# WARNING: do not allow arbitrary user data for date cutoff and lineage set, SQL injection risk

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT organism_taxa.reference_accession, variant_sites.ref_start, variant_sites.ref_end, \
fasta_records.genbank_accession, lineages.name, fasta_records.date_collected, variant_sites.variant_type, \
variant_sites.variant FROM variant_sites \
INNER JOIN fasta_records ON fasta_records.id=variant_sites.fasta_record_id \
INNER JOIN lineage_calls ON fasta_records.lineage_call_id=lineage_calls.id \
INNER JOIN lineages ON lineage_calls.lineage_id=lineages.id \
INNER JOIN organism_taxa ON organism_taxa.id=variant_sites.organism_taxon_id \
INNER JOIN organisms ON organisms.id=organism_taxa.organism_id \
WHERE (date_collected >= '$1' \
OR (date_collected IS NULL AND date_submitted >= '$1')) AND organisms.slug='$2' \
ORDER BY variant_sites.ref_start, variant_sites.ref_end, variant_sites.variant_type, variant_sites.variant;" --csv -t \
| tr "," "\t" | grep -v -E "[[:blank:]][0-9]+N"
