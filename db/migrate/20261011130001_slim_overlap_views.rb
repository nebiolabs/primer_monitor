# frozen_string_literal: true

# The nightly refresh materialized every variant ever seen in variant_overlaps, then cut it to the last 168 days by
# inner-joining time_counts in oligo_variant_overlaps. Apply that window in variant_overlaps instead, so the refresh
# only touches recent sequences, and drop what nothing reads any more: initial_score, counts, time_counts and
# date_counts (only a retired Tableau workbook used them). Results of oligo_variant_overlaps and
# identify_primers_for_notifications are unchanged apart from the dropped count columns.
#
# The views come back empty (WITH NO DATA) until the next nightly refresh; run
# lib/backend_scripts/primer-monitor.sh's REFRESH list by hand after deploying to fill them sooner.
class SlimOverlapViews < ActiveRecord::Migration[8.1] # rubocop:disable Metrics/ClassLength -- view SQL inline
  def up # rubocop:disable Metrics/MethodLength
    execute <<~SQL.squish
      DROP MATERIALIZED VIEW IF EXISTS identify_primers_for_notifications;
      DROP MATERIALIZED VIEW IF EXISTS initial_score;
      DROP MATERIALIZED VIEW IF EXISTS oligo_variant_overlaps;
      DROP MATERIALIZED VIEW IF EXISTS variant_overlaps;
      DROP MATERIALIZED VIEW IF EXISTS time_counts;
      DROP MATERIALIZED VIEW IF EXISTS counts;
      DROP VIEW IF EXISTS date_counts;
    SQL

    create_view 'variant_overlaps', <<~SQL.squish
      SELECT oligos.id AS oligo_id,
                oligos.name AS oligo_name,
                oligo_alignment_positions.ref_start AS oligo_start,
                oligo_alignment_positions.ref_end AS oligo_end,
                oligo_alignment_positions.organism_taxon_id AS alignment_organism_taxon_id,
                oligos.short_name AS oligo_short_name,
                oligos.locus AS oligo_locus_name,
                oligos.category AS oligo_primer_type,
                primer_sets.id AS primer_set_id,
                primer_sets.name AS primer_set_name,
                variant_sites.id AS variant_id,
                variant_sites.variant_type,
                variant_sites.variant,
                variant_sites.ref_start AS variant_start,
                variant_sites.ref_end AS variant_end,
                COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
                COALESCE(detailed_geo_locations.subregion, ''::character varying) AS subregion,
                COALESCE(detailed_geo_locations.division, ''::character varying) AS division,
                COALESCE(detailed_geo_locations.subdivision, ''::character varying) AS subdivision,
                detailed_geo_locations.id AS detailed_geo_location_id,
                COALESCE(fasta_records.date_collected, '1900-01-01'::date) AS date_collected
         FROM variant_sites
                  JOIN fasta_records ON variant_sites.fasta_record_id = fasta_records.id
                  JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
                  JOIN oligo_alignment_positions ON NOT (oligo_alignment_positions.ref_start >= variant_sites.ref_end OR oligo_alignment_positions.ref_end <= variant_sites.ref_start) AND fasta_records.organism_taxon_id = oligo_alignment_positions.organism_taxon_id
                  JOIN oligos ON oligos.id = oligo_alignment_positions.oligo_id
                  JOIN primer_sets ON oligos.primer_set_id = primer_sets.id
         WHERE variant_sites.usable_del_or_snp = true AND fasta_records.date_collected > (SELECT max(recent.date_collected) - '168 days'::interval FROM fasta_records recent)
         UNION ALL
         SELECT oligos.id AS oligo_id,
                oligos.name AS oligo_name,
                oligo_alignment_positions.ref_start AS oligo_start,
                oligo_alignment_positions.ref_end AS oligo_end,
                oligo_alignment_positions.organism_taxon_id AS alignment_organism_taxon_id,
                oligos.short_name AS oligo_short_name,
                oligos.locus AS oligo_locus_name,
                oligos.category AS oligo_primer_type,
                primer_sets.id AS primer_set_id,
                primer_sets.name AS primer_set_name,
                variant_sites.id AS variant_id,
                variant_sites.variant_type,
                variant_sites.variant,
                variant_sites.ref_start AS variant_start,
                variant_sites.ref_end AS variant_end,
                COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
                COALESCE(detailed_geo_locations.subregion, ''::character varying) AS subregion,
                COALESCE(detailed_geo_locations.division, ''::character varying) AS division,
                COALESCE(detailed_geo_locations.subdivision, ''::character varying) AS subdivision,
                detailed_geo_locations.id AS detailed_geo_location_id,
                COALESCE(fasta_records.date_collected, '1900-01-01'::date) AS date_collected
         FROM variant_sites
                  JOIN fasta_records ON variant_sites.fasta_record_id = fasta_records.id
                  JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
                  JOIN oligo_alignment_positions ON NOT (oligo_alignment_positions.ref_start >= variant_sites.ref_end OR oligo_alignment_positions.ref_end <= variant_sites.ref_start) AND fasta_records.organism_taxon_id = oligo_alignment_positions.organism_taxon_id
                  JOIN oligos ON oligos.id = oligo_alignment_positions.oligo_id
                  JOIN primer_sets ON oligos.primer_set_id = primer_sets.id
         WHERE variant_sites.usable_insertion = true AND fasta_records.date_collected > (SELECT max(recent.date_collected) - '168 days'::interval FROM fasta_records recent)
    SQL
    create_view 'oligo_variant_overlaps', <<~SQL.squish
      SELECT variant_overlaps.oligo_id,
             variant_overlaps.oligo_name,
             variant_overlaps.oligo_start,
             variant_overlaps.oligo_end,
             variant_overlaps.oligo_short_name,
             variant_overlaps.oligo_locus_name,
             variant_overlaps.oligo_primer_type,
             variant_overlaps.primer_set_id,
             variant_overlaps.primer_set_name,
             variant_overlaps.variant_id,
             variant_overlaps.variant_type,
             variant_overlaps.variant,
             variant_overlaps.variant_start,
             variant_overlaps.variant_end,
             variant_overlaps.region,
             variant_overlaps.subregion,
             variant_overlaps.division,
             variant_overlaps.subdivision,
             variant_overlaps.detailed_geo_location_id,
             variant_overlaps.date_collected,
             generate_series(lower((numrange((variant_overlaps.oligo_start)::numeric, (variant_overlaps.oligo_end)::numeric) *
                                    numrange((variant_overlaps.variant_start)::numeric,
                                             (variant_overlaps.variant_end)::numeric))),
                             (upper((numrange((variant_overlaps.oligo_start)::numeric,
                                              (variant_overlaps.oligo_end)::numeric) *
                                     numrange((variant_overlaps.variant_start)::numeric,
                                              (variant_overlaps.variant_end)::numeric))) - (1)::numeric)) AS coords
      FROM public.variant_overlaps variant_overlaps
      WHERE ((variant_overlaps.variant_type)::text <> 'I'::text)
      UNION ALL
      SELECT variant_overlaps.oligo_id,
             variant_overlaps.oligo_name,
             variant_overlaps.oligo_start,
             variant_overlaps.oligo_end,
             variant_overlaps.oligo_short_name,
             variant_overlaps.oligo_locus_name,
             variant_overlaps.oligo_primer_type,
             variant_overlaps.primer_set_id,
             variant_overlaps.primer_set_name,
             variant_overlaps.variant_id,
             variant_overlaps.variant_type,
             variant_overlaps.variant,
             variant_overlaps.variant_start,
             variant_overlaps.variant_end,
             variant_overlaps.region,
             variant_overlaps.subregion,
             variant_overlaps.division,
             variant_overlaps.subdivision,
             variant_overlaps.detailed_geo_location_id,
             variant_overlaps.date_collected,
             variant_overlaps.variant_start AS coords
      FROM public.variant_overlaps
      WHERE ((variant_overlaps.variant_type)::text = 'I'::text)
    SQL
    %w[oligo_id date_collected detailed_geo_location_id].each do |column|
      execute "CREATE INDEX ON oligo_variant_overlaps (#{column})"
    end
    create_view 'identify_primers_for_notifications', <<~SQL.squish
      WITH subscribed_locations AS (
               SELECT DISTINCT join_subscribed_location_to_ids.user_id,
                  join_subscribed_location_to_ids.detailed_geo_location_id
                 FROM public.join_subscribed_location_to_ids
              ), first_query AS (
               SELECT primer_set_subscriptions.user_id,
                  primer_set_subscriptions.primer_set_id,
                  subscribed_locations.detailed_geo_location_id,
                  primer_sets.name AS set_name,
                  oligos.name AS primer_name,
                  oligos.id AS oligo_id,
                  users.lookback_days,
                  users.variant_fraction_threshold,
                  oligo_variant_overlaps.region,
                  oligo_variant_overlaps.subregion,
                  oligo_variant_overlaps.division,
                  oligo_variant_overlaps.subdivision,
                  oligo_variant_overlaps.detailed_geo_location_id AS unused_id,
                  oligo_variant_overlaps.coords,
                  count(oligo_variant_overlaps.variant_id) AS variant_count
                 FROM (((((public.primer_set_subscriptions
                   JOIN public.primer_sets ON ((primer_sets.id = primer_set_subscriptions.primer_set_id)))
                   JOIN public.oligos ON ((primer_sets.id = oligos.primer_set_id)))
                   JOIN public.oligo_variant_overlaps ON ((oligo_variant_overlaps.oligo_id = oligos.id)))
                   JOIN subscribed_locations ON (((subscribed_locations.user_id = primer_set_subscriptions.user_id) AND (subscribed_locations.detailed_geo_location_id = oligo_variant_overlaps.detailed_geo_location_id))))
                   JOIN public.users ON ((users.id = primer_set_subscriptions.user_id)))
                WHERE ((oligo_variant_overlaps.date_collected >= (CURRENT_DATE - users.lookback_days)) AND (primer_set_subscriptions.active = true))
                GROUP BY primer_set_subscriptions.user_id, primer_set_subscriptions.primer_set_id, primer_sets.name, oligos.id, oligos.name, subscribed_locations.detailed_geo_location_id, users.lookback_days, users.variant_fraction_threshold, oligo_variant_overlaps.region, oligo_variant_overlaps.subregion, oligo_variant_overlaps.division, oligo_variant_overlaps.subdivision, oligo_variant_overlaps.coords, oligo_variant_overlaps.detailed_geo_location_id
              ), total_sequences_for_denominator AS (
               SELECT fasta_records.detailed_geo_location_id,
                  count(fasta_records.id) AS records_count,
                  location_lookbacks.lookback_days
                 FROM (public.fasta_records
                   JOIN ( SELECT DISTINCT subscribed_locations.detailed_geo_location_id,
                          users.lookback_days
                         FROM (subscribed_locations
                           JOIN public.users ON ((users.id = subscribed_locations.user_id)))) location_lookbacks ON ((location_lookbacks.detailed_geo_location_id = fasta_records.detailed_geo_location_id)))
                WHERE (fasta_records.date_collected >= (CURRENT_DATE - location_lookbacks.lookback_days))
                GROUP BY fasta_records.detailed_geo_location_id, location_lookbacks.lookback_days
               HAVING (count(fasta_records.id) >= 20)
              )
       SELECT first_query.user_id,
          first_query.primer_set_id,
          first_query.set_name,
          first_query.oligo_id,
          first_query.primer_name,
          first_query.region,
          first_query.subregion,
          first_query.division,
          first_query.subdivision,
          first_query.coords,
          first_query.variant_count,
          total_sequences_for_denominator.records_count,
          first_query.variant_fraction_threshold,
          total_sequences_for_denominator.detailed_geo_location_id,
          ((first_query.variant_count)::numeric / (total_sequences_for_denominator.records_count)::numeric) AS fraction_variant
         FROM (first_query
           JOIN total_sequences_for_denominator ON (((total_sequences_for_denominator.detailed_geo_location_id = first_query.detailed_geo_location_id) AND (total_sequences_for_denominator.lookback_days = first_query.lookback_days))))
        WHERE ((((first_query.variant_count)::numeric / (total_sequences_for_denominator.records_count)::numeric))::double precision >= first_query.variant_fraction_threshold)
    SQL
  end

  def down
    # the previous definitions are in db/views and db/structure.sql at the parent commit
    raise ActiveRecord::IrreversibleMigration
  end

  private

  def create_view(name, definition)
    execute "CREATE MATERIALIZED VIEW #{name} AS #{definition} WITH NO DATA"
    execute "GRANT SELECT ON #{name} TO primer_monitor_ro"
  end
end
