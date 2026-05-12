#!/usr/bin/env bash

# Stores variant-primer overlaps for one lineage group into the DB.
#
# Usage: store_variant_overlaps.sh <organism_slug> <lineage_group_key> <variants_bed_path>
#
# Requires: DB_HOST, DB_NAME, DB_USER (write-capable), DB_PASSWORD env vars exported by caller.
# variants_bed_path: the lineage_variants/{lineage}.bed file produced by count_variants.sh
#   Columns: chrom, ref_start, ref_end, variant (allele only), frequency_pct
# variant_type is looked up from variant_sites rather than parsed from the filename.
#
# WARNING: organism_slug, lineage_group_key, and variants_bed_path are interpolated directly
# into SQL. Do not pass arbitrary user input to this script.

set -e

organism_slug="$1"
lineage_group_key="$2"
variants_bed="$3"

if [[ -z "$organism_slug" || -z "$lineage_group_key" || -z "$variants_bed" ]]; then
  echo "Usage: store_variant_overlaps.sh <organism_slug> <lineage_group_key> <variants_bed_path>" >&2
  exit 1
fi

sql_file=$(mktemp /tmp/store_overlaps_XXXXXX.sql)
trap "rm -f '$sql_file'" EXIT

cat > "$sql_file" <<SQL
\set ON_ERROR_STOP on

CREATE TEMP TABLE tmp_variants (
  chrom        varchar,
  ref_start    integer,
  ref_end      integer,
  variant_name varchar,
  frequency_pct float
) ON COMMIT DROP;

\COPY tmp_variants FROM '${variants_bed}' WITH (FORMAT csv, DELIMITER E'\\t', HEADER false)

-- Remove stale rows from previous pipeline runs for this organism + lineage group
DELETE FROM lineage_variant_primer_overlaps
WHERE organism_id = (SELECT id FROM organisms WHERE slug = '${organism_slug}')
  AND lineage_group_key = '${lineage_group_key}';

-- Compute first/last seen globally for each distinct variant position.
CREATE TEMP TABLE tmp_seen AS
SELECT
  vs.ref_start,
  vs.ref_end,
  vs.variant_type,
  vs.variant,
  MIN(COALESCE(fr.date_collected, fr.date_submitted))
    FILTER (WHERE COALESCE(fr.date_collected, fr.date_submitted) IS NOT NULL) AS first_date,
  MAX(COALESCE(fr.date_collected, fr.date_submitted))
    FILTER (WHERE COALESCE(fr.date_collected, fr.date_submitted) IS NOT NULL) AS last_date
FROM tmp_variants tv
JOIN organism_taxa ot2 ON ot2.reference_accession = tv.chrom
JOIN variant_sites vs  ON vs.ref_start = tv.ref_start
                      AND vs.ref_end   = tv.ref_end
                      AND vs.variant   = tv.variant_name
                      AND vs.organism_taxon_id = ot2.id
JOIN fasta_records fr  ON fr.id = vs.fasta_record_id
GROUP BY vs.ref_start, vs.ref_end, vs.variant_type, vs.variant;

-- Enrich first_date with lineage + location.
CREATE TEMP TABLE tmp_first AS
SELECT DISTINCT ON (vs.ref_start, vs.variant_type, vs.variant)
  vs.ref_start, vs.variant_type, vs.variant,
  l.name AS lineage,
  COALESCE(dga.region, '') || CASE WHEN dga.division IS NOT NULL THEN ' / ' || dga.division ELSE '' END AS location
FROM tmp_seen s
JOIN variant_sites vs  ON vs.ref_start = s.ref_start AND vs.ref_end = s.ref_end
                      AND vs.variant_type = s.variant_type AND vs.variant = s.variant
JOIN fasta_records fr  ON fr.id = vs.fasta_record_id
                      AND COALESCE(fr.date_collected, fr.date_submitted) = s.first_date
JOIN lineage_calls lc  ON lc.id = fr.lineage_call_id
JOIN lineages l        ON l.id  = lc.lineage_id
JOIN detailed_geo_locations dgl   ON dgl.id = fr.detailed_geo_location_id
JOIN detailed_geo_location_aliases dga ON dga.id = dgl.detailed_geo_location_alias_id
ORDER BY vs.ref_start, vs.variant_type, vs.variant, fr.id ASC;

-- Enrich last_date with lineage + location.
CREATE TEMP TABLE tmp_last AS
SELECT DISTINCT ON (vs.ref_start, vs.variant_type, vs.variant)
  vs.ref_start, vs.variant_type, vs.variant,
  l.name AS lineage,
  COALESCE(dga.region, '') || CASE WHEN dga.division IS NOT NULL THEN ' / ' || dga.division ELSE '' END AS location
FROM tmp_seen s
JOIN variant_sites vs  ON vs.ref_start = s.ref_start AND vs.ref_end = s.ref_end
                      AND vs.variant_type = s.variant_type AND vs.variant = s.variant
JOIN fasta_records fr  ON fr.id = vs.fasta_record_id
                      AND COALESCE(fr.date_collected, fr.date_submitted) = s.last_date
JOIN lineage_calls lc  ON lc.id = fr.lineage_call_id
JOIN lineages l        ON l.id  = lc.lineage_id
JOIN detailed_geo_locations dgl   ON dgl.id = fr.detailed_geo_location_id
JOIN detailed_geo_location_aliases dga ON dga.id = dgl.detailed_geo_location_alias_id
ORDER BY vs.ref_start, vs.variant_type, vs.variant, fr.id DESC;

-- Insert overlaps: variant positions x oligos aligned to the same reference.
-- variant_type is looked up from variant_sites; first/last seen pre-computed above.
INSERT INTO lineage_variant_primer_overlaps
  (organism_id, lineage_group_key, ref_start, ref_end,
   variant_type, variant, frequency_pct, oligo_id,
   first_seen_date, first_seen_lineage, first_seen_location,
   last_seen_date,  last_seen_lineage,  last_seen_location)
SELECT DISTINCT
  org.id,
  '${lineage_group_key}',
  tv.ref_start,
  tv.ref_end,
  vs.variant_type,
  tv.variant_name,
  tv.frequency_pct,
  o.id,
  s.first_date,  tf.lineage,  tf.location,
  s.last_date,   tl.lineage,  tl.location
FROM tmp_variants tv
JOIN organism_taxa ot  ON ot.reference_accession = tv.chrom
JOIN organisms org     ON org.id = ot.organism_id AND org.slug = '${organism_slug}'
JOIN variant_sites vs  ON vs.ref_start = tv.ref_start
                      AND vs.ref_end   = tv.ref_end
                      AND vs.variant   = tv.variant_name
                      AND vs.organism_taxon_id = ot.id
JOIN tmp_seen s        ON s.ref_start = vs.ref_start AND s.variant_type = vs.variant_type AND s.variant = vs.variant
LEFT JOIN tmp_first tf ON tf.ref_start = vs.ref_start AND tf.variant_type = vs.variant_type AND tf.variant = vs.variant
LEFT JOIN tmp_last  tl ON tl.ref_start = vs.ref_start AND tl.variant_type = vs.variant_type AND tl.variant = vs.variant
JOIN oligo_alignment_positions oap
  ON  oap.organism_taxon_id = ot.id
  AND NOT (oap.ref_start >= tv.ref_end OR oap.ref_end <= tv.ref_start)
JOIN oligos o          ON o.id = oap.oligo_id
JOIN primer_sets ps    ON ps.id = o.primer_set_id AND ps.organism_id = org.id
ON CONFLICT DO NOTHING;
SQL

PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" -f "$sql_file"
