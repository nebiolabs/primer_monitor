#!/usr/bin/env bash

# Stores variant-primer overlaps for one lineage group into the DB.
#
# Usage: store_variant_overlaps.sh <organism_slug> <lineage_group_key> <variants_bed_path>
#
# Requires: DB_HOST, DB_NAME, DB_USER (write-capable) env vars exported by caller.
# variants_bed_path: the lineage_variants/{lineage}.bed file produced by count_variants.sh
#   Columns: chrom, ref_start, ref_end, variant_name (type/allele), frequency_pct

set -e

organism_slug="$1"
lineage_group_key="$2"
variants_bed="$3"

if [[ -z "$organism_slug" || -z "$lineage_group_key" || -z "$variants_bed" ]]; then
  echo "Usage: store_variant_overlaps.sh <organism_slug> <lineage_group_key> <variants_bed_path>" >&2
  exit 1
fi

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" \
  -v "organism_slug=$organism_slug" \
  -v "lineage_key=$lineage_group_key" \
  -v "variants_bed=$variants_bed" \
  -f - <<'PSQL'

CREATE TEMP TABLE tmp_variants (
  chrom        varchar,
  ref_start    integer,
  ref_end      integer,
  variant_name varchar,
  frequency_pct float
) ON COMMIT DROP;

\COPY tmp_variants FROM :'variants_bed' WITH (FORMAT csv, DELIMITER E'\t', HEADER false)

-- Remove stale rows from previous pipeline runs for this organism + lineage group
DELETE FROM lineage_variant_primer_overlaps
WHERE organism_id = (SELECT id FROM organisms WHERE slug = :'organism_slug')
  AND lineage_group_key = :'lineage_key';

-- Insert overlaps: variant positions x oligos aligned to the same reference
INSERT INTO lineage_variant_primer_overlaps
  (organism_id, lineage_group_key, ref_start, ref_end,
   variant_type, variant, frequency_pct, oligo_id)
SELECT
  org.id,
  :'lineage_key',
  tv.ref_start,
  tv.ref_end,
  split_part(tv.variant_name, '/', 1),
  split_part(tv.variant_name, '/', 2),
  tv.frequency_pct,
  o.id
FROM tmp_variants tv
JOIN organism_taxa ot  ON ot.reference_accession = tv.chrom
JOIN organisms org     ON org.id = ot.organism_id AND org.slug = :'organism_slug'
JOIN oligo_alignment_positions oap
  ON  oap.organism_taxon_id = ot.id
  AND NOT (oap.ref_start >= tv.ref_end OR oap.ref_end <= tv.ref_start)
JOIN oligos o          ON o.id = oap.oligo_id
JOIN primer_sets ps    ON ps.id = o.primer_set_id AND ps.organism_id = org.id
ON CONFLICT DO NOTHING;

PSQL
