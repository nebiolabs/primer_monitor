#!/usr/bin/env bash

if (($# < 2)); then
    echo "usage: extract_all_variants.sh <date cutoff> <threads>" >&2;
    exit 1;
fi
# you need to export DB_HOST, DB_NAME, and DB_USER before running this

# WARNING: do not allow arbitrary user data for date cutoff and lineage set, SQL injection risk

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" -c "SELECT variant_sites.ref_start, variant_sites.ref_end, fasta_records.strain, fasta_records.pangolin_lineage, fasta_records.date_collected, variant_sites.variant_type, variant_sites.variant FROM variant_sites INNER JOIN fasta_records ON fasta_records.id=variant_sites.fasta_record_id WHERE (date_collected >= '$1' OR (date_collected IS NULL AND date_submitted >= '$1'));" --csv -t | tr "," "\t" | sed -E "s/^/NC_045512.2\t/g" | sort --parallel="$2" -k2n -k3n