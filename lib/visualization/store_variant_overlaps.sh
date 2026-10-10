#!/usr/bin/env bash

# Stores the variant-primer overlaps of all of an organism's lineage groups into the DB, replacing the old ones.
#
# Usage: store_variant_overlaps.sh <organism_slug> <lineage_variants_dir> <lineage_group_key>...
#
# Requires: DB_HOST, DB_NAME, DB_USER (write-capable), DB_PASSWORD env vars exported by caller.
# <lineage_variants_dir>/<lineage_group_key>.bed: produced by count_variants.sh
#   Columns: chrom, ref_start, ref_end, variant (allele only), frequency_pct
#
# First/last seen are over all sequences, whatever their lineage, so they are computed once per variant rather than
# once per lineage group: variant_sites holds hundreds of millions of rows, and a common variant matches millions.
# Only variants overlapping an aligned primer are looked up, since only those are stored.
#
# WARNING: organism_slug, lineage group keys and paths are interpolated directly
# into SQL. Do not pass arbitrary user input to this script.

set -e

organism_slug="$1"
variants_dir="$2"
shift 2

if [[ -z "$organism_slug" || -z "$variants_dir" || $# -eq 0 ]]; then
  echo "Usage: store_variant_overlaps.sh <organism_slug> <lineage_variants_dir> <lineage_group_key>..." >&2
  exit 1
fi

sql_file=$(mktemp /tmp/store_overlaps_XXXXXX.sql)
trap "rm -f '$sql_file'" EXIT

copy_commands=""
for lineage_group_key in "$@"; do
  copy_commands+="\\COPY tmp_lineage_variants FROM '${variants_dir}/${lineage_group_key}.bed' WITH (FORMAT csv, DELIMITER E'\\t', HEADER false)
INSERT INTO tmp_variants SELECT *, '${lineage_group_key}' FROM tmp_lineage_variants;
TRUNCATE tmp_lineage_variants;
"
done

cat > "$sql_file" <<SQL
\set ON_ERROR_STOP on

CREATE TEMP TABLE tmp_lineage_variants (
  chrom         varchar,
  ref_start     integer,
  ref_end       integer,
  variant_name  varchar,
  frequency_pct float
);
CREATE TEMP TABLE tmp_variants (LIKE tmp_lineage_variants, lineage_group_key varchar);

${copy_commands}
ANALYZE tmp_variants;

-- the distinct variants that overlap one of the organism's aligned primers
CREATE TEMP TABLE tmp_overlapping AS
SELECT DISTINCT ot.id AS organism_taxon_id, tv.ref_start, tv.ref_end, tv.variant_name AS variant
FROM tmp_variants tv
JOIN organism_taxa ot  ON ot.reference_accession = tv.chrom
JOIN organisms org     ON org.id = ot.organism_id AND org.slug = '${organism_slug}'
WHERE EXISTS (
  SELECT 1 FROM oligo_alignment_positions oap
  JOIN oligos o       ON o.id = oap.oligo_id
  JOIN primer_sets ps ON ps.id = o.primer_set_id AND ps.organism_id = org.id
  WHERE oap.organism_taxon_id = ot.id AND NOT (oap.ref_start >= tv.ref_end OR oap.ref_end <= tv.ref_start)
);
ANALYZE tmp_overlapping;

-- One pass over each variant's sequences. A key of date then id compares like (date, id), so min/max pick the
-- earliest sequence (lowest id on ties) and the latest (highest id) without sorting millions of rows.
CREATE TEMP TABLE tmp_seen AS
SELECT v.organism_taxon_id, v.ref_start, v.ref_end, vs.variant_type, v.variant, max(vs.ref) AS ref,
  min(to_char(COALESCE(fr.date_collected, fr.date_submitted), 'YYYYMMDD') || lpad(fr.id::text, 19, '0')) AS first_key,
  max(to_char(COALESCE(fr.date_collected, fr.date_submitted), 'YYYYMMDD') || lpad(fr.id::text, 19, '0')) AS last_key
FROM tmp_overlapping v
JOIN variant_sites vs  ON vs.organism_taxon_id = v.organism_taxon_id
                      AND vs.ref_start = v.ref_start AND vs.ref_end = v.ref_end AND vs.variant = v.variant
JOIN fasta_records fr  ON fr.id = vs.fasta_record_id
GROUP BY v.organism_taxon_id, v.ref_start, v.ref_end, vs.variant_type, v.variant;

-- the lineage and location of the first and last sequences
CREATE TEMP TABLE tmp_seen_info AS
SELECT s.organism_taxon_id, s.ref_start, s.ref_end, s.variant_type, s.variant, s.ref,
  to_date(left(s.first_key, 8), 'YYYYMMDD') AS first_date, fl.name AS first_lineage,
  COALESCE(fa.region, '') || CASE WHEN fa.division IS NOT NULL THEN ' / ' || fa.division ELSE '' END AS first_location,
  to_date(left(s.last_key, 8), 'YYYYMMDD') AS last_date, ll.name AS last_lineage,
  COALESCE(la.region, '') || CASE WHEN la.division IS NOT NULL THEN ' / ' || la.division ELSE '' END AS last_location
FROM tmp_seen s
LEFT JOIN fasta_records ff                ON ff.id = substr(s.first_key, 9)::bigint
LEFT JOIN lineage_calls fc                ON fc.id = ff.lineage_call_id
LEFT JOIN lineages fl                     ON fl.id = fc.lineage_id
LEFT JOIN detailed_geo_locations fg       ON fg.id = ff.detailed_geo_location_id
LEFT JOIN detailed_geo_location_aliases fa ON fa.id = fg.detailed_geo_location_alias_id
LEFT JOIN fasta_records lf                ON lf.id = substr(s.last_key, 9)::bigint
LEFT JOIN lineage_calls lc                ON lc.id = lf.lineage_call_id
LEFT JOIN lineages ll                     ON ll.id = lc.lineage_id
LEFT JOIN detailed_geo_locations lg       ON lg.id = lf.detailed_geo_location_id
LEFT JOIN detailed_geo_location_aliases la ON la.id = lg.detailed_geo_location_alias_id;

-- replace all of the organism's overlaps, so lineage groups no longer shown don't linger
DELETE FROM lineage_variant_primer_overlaps
WHERE organism_id = (SELECT id FROM organisms WHERE slug = '${organism_slug}');

INSERT INTO lineage_variant_primer_overlaps
  (organism_id, lineage_group_key, ref_start, ref_end,
   variant_type, variant, ref, frequency_pct, oligo_id,
   first_seen_date, first_seen_lineage, first_seen_location,
   last_seen_date,  last_seen_lineage,  last_seen_location)
SELECT DISTINCT
  org.id, tv.lineage_group_key, tv.ref_start, tv.ref_end,
  s.variant_type, tv.variant_name, s.ref, tv.frequency_pct, o.id,
  s.first_date, s.first_lineage, s.first_location,
  s.last_date,  s.last_lineage,  s.last_location
FROM tmp_variants tv
JOIN organism_taxa ot  ON ot.reference_accession = tv.chrom
JOIN organisms org     ON org.id = ot.organism_id AND org.slug = '${organism_slug}'
JOIN tmp_seen_info s   ON s.organism_taxon_id = ot.id AND s.ref_start = tv.ref_start
                      AND s.ref_end = tv.ref_end AND s.variant = tv.variant_name
JOIN oligo_alignment_positions oap
  ON  oap.organism_taxon_id = ot.id
  AND NOT (oap.ref_start >= tv.ref_end OR oap.ref_end <= tv.ref_start)
JOIN oligos o          ON o.id = oap.oligo_id
JOIN primer_sets ps    ON ps.id = o.primer_set_id AND ps.organism_id = org.id
ON CONFLICT DO NOTHING;
SQL

# one transaction: the page never sees the organism's overlaps half replaced
PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" --single-transaction -f "$sql_file"
