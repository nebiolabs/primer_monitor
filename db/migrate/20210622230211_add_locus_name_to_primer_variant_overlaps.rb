class AddLocusNameToPrimerVariantOverlaps < ActiveRecord::Migration[6.1]
  def up
    execute <<~SQL
      DROP MATERIALIZED VIEW IF EXISTS identify_primers_for_notifications;
      DROP MATERIALIZED VIEW IF EXISTS initial_score;
      DROP MATERIALIZED VIEW IF EXISTS oligo_variant_overlaps;    
      DROP MATERIALIZED VIEW IF EXISTS variant_overlaps;  

      CREATE MATERIALIZED VIEW IF NOT EXISTS variant_overlaps AS
          SELECT oligos.id AS oligo_id,
                oligos.name AS oligo_name,
                oligos.ref_start AS oligo_start,
                oligos.ref_end AS oligo_end,
                oligos.short_name as oligo_short_name,
                oligos.locus as oligo_locus_name,
                oligos.category as oligo_primer_type,
                primer_sets.id as primer_set_id,
                primer_sets.name as primer_set_name,
                variant_sites.id AS variant_id,
                variant_sites.variant_type,
                variant_sites.variant,
                variant_sites.ref_start AS variant_start,
                variant_sites.ref_end AS variant_end,
                coalesce(detailed_geo_locations.region, '') as region,
                coalesce(detailed_geo_locations.subregion, '') as subregion,
                coalesce(detailed_geo_locations.division, '') as division,
                coalesce(detailed_geo_locations.subdivision, '') as subdivision,
                detailed_geo_locations.id AS detailed_geo_location_id,
                coalesce(fasta_records.date_collected, '1900-01-01') as date_collected
              FROM variant_sites
              JOIN fasta_records ON variant_sites.fasta_record_id = fasta_records.id
              JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
              JOIN oligos on NOT (oligos.ref_start >= variant_sites.ref_end or oligos.ref_end <= variant_sites.ref_start )
              JOIN primer_sets on oligos.primer_set_id = primer_sets.id
        WHERE usable_del_or_snp = true 
        UNION ALL 
        SELECT oligos.id AS oligo_id,
              oligos.name AS oligo_name,
              oligos.ref_start AS oligo_start,
              oligos.ref_end AS oligo_end,
              oligos.short_name as oligo_short_name,
              oligos.locus as oligo_locus_name,
              oligos.category as oligo_primer_type,
              primer_sets.id as primer_set_id,
              primer_sets.name as primer_set_name,
              variant_sites.id AS variant_id,
              variant_sites.variant_type,
              variant_sites.variant,
              variant_sites.ref_start AS variant_start,
              variant_sites.ref_end AS variant_end,
              coalesce(detailed_geo_locations.region, '') as region,
              coalesce(detailed_geo_locations.subregion, '') as subregion,
              coalesce(detailed_geo_locations.division, '') as division,
              coalesce(detailed_geo_locations.subdivision, '') as subdivision,
              detailed_geo_locations.id AS detailed_geo_location_id,
              coalesce(fasta_records.date_collected, '1900-01-01') as date_collected
            FROM variant_sites
            JOIN fasta_records ON variant_sites.fasta_record_id = fasta_records.id
            JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
            JOIN oligos ON variant_sites.ref_start > oligos.ref_start AND variant_sites.ref_start <= oligos.ref_end 
            JOIN primer_sets on oligos.primer_set_id = primer_sets.id
        WHERE usable_insertion = true 
      WITH DATA;

      CREATE INDEX ON variant_overlaps(region, subregion, division, subdivision, date_collected);

      CREATE MATERIALIZED VIEW IF NOT EXISTS oligo_variant_overlaps AS
        SELECT variant_overlaps.*, counts.region_count, counts.region_subregion_count, counts.region_subregion_division_count, counts.region_subregion_division_subdivision_count,
          time_counts.region_time_count, time_counts.region_subregion_time_count, time_counts.region_subregion_division_time_count, time_counts.region_subregion_division_subdivision_time_count,
          GENERATE_SERIES(
          LOWER(numrange(variant_overlaps.oligo_start, variant_overlaps.oligo_end) * numrange(variant_overlaps.variant_start, variant_overlaps.variant_end)),
          UPPER(numrange(variant_overlaps.oligo_start, variant_overlaps.oligo_end) * numrange(variant_overlaps.variant_start, variant_overlaps.variant_end)) - 1)
        AS coords FROM variant_overlaps AS variant_overlaps
        INNER JOIN counts on counts.region = variant_overlaps.region
          AND counts.subregion = variant_overlaps.subregion
          AND counts.division = variant_overlaps.division
          AND counts.subdivision = variant_overlaps.subdivision
        INNER JOIN time_counts on time_counts.region = variant_overlaps.region
          AND time_counts.subregion = variant_overlaps.subregion
          AND time_counts.division = variant_overlaps.division
          AND time_counts.subdivision = variant_overlaps.subdivision
          AND time_counts.date_collected = variant_overlaps.date_collected
        WHERE variant_type <> 'I'
        UNION ALL
        SELECT variant_overlaps.*, counts.region_count, counts.region_subregion_count, counts.region_subregion_division_count, counts.region_subregion_division_subdivision_count,
          time_counts.region_time_count, time_counts.region_subregion_time_count, time_counts.region_subregion_division_time_count, time_counts.region_subregion_division_subdivision_time_count,
          variant_overlaps.variant_start as coords
        FROM variant_overlaps
        INNER JOIN counts on counts.region = variant_overlaps.region
          AND counts.subregion = variant_overlaps.subregion
          AND counts.division = variant_overlaps.division
          AND counts.subdivision = variant_overlaps.subdivision
        INNER JOIN time_counts on time_counts.region = variant_overlaps.region
          AND time_counts.subregion = variant_overlaps.subregion
          AND time_counts.division = variant_overlaps.division
          AND time_counts.subdivision = variant_overlaps.subdivision
          AND time_counts.date_collected = variant_overlaps.date_collected
        WHERE variant_type = 'I'
      WITH DATA;

      CREATE INDEX ON oligo_variant_overlaps(primer_set_name);
      CREATE INDEX ON oligo_variant_overlaps(region, subregion, division, subdivision, date_collected);


     CREATE MATERIALIZED VIEW IF NOT EXISTS identify_primers_for_notifications AS
      WITH first_query AS (
         SELECT primer_set_subscriptions.user_id,
            primer_set_subscriptions.primer_set_id,
            join_subscribed_location_to_ids.detailed_geo_location_id,
            primer_sets.name AS set_name,
            oligos.name AS primer_name,
            oligos.id as oligo_id,
            users.lookback_days,
            users.variant_fraction_threshold,
            oligo_variant_overlaps.region,
            oligo_variant_overlaps.subregion,
            oligo_variant_overlaps.division,
            oligo_variant_overlaps.subdivision,
            oligo_variant_overlaps.detailed_geo_location_id AS unused_id,
            oligo_variant_overlaps.coords,
            count(oligo_variant_overlaps.variant_id) AS variant_count
           FROM primer_set_subscriptions
             JOIN primer_sets ON primer_sets.id = primer_set_subscriptions.primer_set_id
             JOIN oligos ON primer_sets.id = oligos.primer_set_id
             JOIN oligo_variant_overlaps ON oligo_variant_overlaps.oligo_id = oligos.id
             JOIN join_subscribed_location_to_ids ON join_subscribed_location_to_ids.user_id = primer_set_subscriptions.user_id and join_subscribed_location_to_ids.detailed_geo_location_id = oligo_variant_overlaps.detailed_geo_location_id
             JOIN users ON users.id = primer_set_subscriptions.user_id
          WHERE oligo_variant_overlaps.date_collected >= (CURRENT_DATE - users.lookback_days)
          GROUP BY primer_set_subscriptions.user_id, primer_set_subscriptions.primer_set_id, primer_sets.name, oligos.id, oligos.name,
                   join_subscribed_location_to_ids.detailed_geo_location_id,
                   users.lookback_days, users.variant_fraction_threshold,
                   oligo_variant_overlaps.region, oligo_variant_overlaps.subregion, oligo_variant_overlaps.division, oligo_variant_overlaps.subdivision, oligo_variant_overlaps.coords, oligo_variant_overlaps.detailed_geo_location_id
        ), second_query AS (
         SELECT fasta_records.detailed_geo_location_id,
            count(fasta_records.id) AS records_count,
            users.lookback_days
           FROM fasta_records
             JOIN join_subscribed_location_to_ids ON join_subscribed_location_to_ids.detailed_geo_location_id = fasta_records.detailed_geo_location_id
             JOIN users ON join_subscribed_location_to_ids.user_id = users.id
          WHERE fasta_records.date_collected >= (CURRENT_DATE - users.lookback_days)
          GROUP BY fasta_records.detailed_geo_location_id, users.lookback_days
         HAVING count(fasta_records.id) >= 20
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
            first_query.variant_fraction_threshold,
            second_query.detailed_geo_location_id,
            second_query.records_count,
            first_query.variant_count::numeric / second_query.records_count::numeric AS fraction_variant
        FROM first_query
          JOIN second_query ON second_query.detailed_geo_location_id = first_query.detailed_geo_location_id
        WHERE (first_query.variant_count::numeric / second_query.records_count::numeric)::double precision >= first_query.variant_fraction_threshold
        WITH DATA;

      CREATE MATERIALIZED VIEW IF NOT EXISTS initial_score AS
        WITH all_combos AS (
          SELECT oligos_1.oligo_name,
            oligos_1.locus,
            oligos_1.primer_set_id,
            fa_1.detailed_geo_location_id,
            fa_1.date_collected
            FROM ( SELECT DISTINCT oligos.name AS oligo_name,
                    oligos.locus,
                    oligos.primer_set_id
                    FROM oligos) oligos_1
              CROSS JOIN ( SELECT DISTINCT fasta_records.detailed_geo_location_id,
                    fasta_records.date_collected
                    FROM fasta_records
                  WHERE fasta_records.date_collected >= '2020-01-01'::date) fa_1
        ), variants AS (
          SELECT count(variant_sites.fasta_record_id) AS observed_count,
            ovo.detailed_geo_location_id,
            count(*)::double precision / (5::double precision * count(variant_sites.fasta_record_id)::double precision) AS three_p_score,
            oligos.locus,
            ovo.primer_set_id,
            ovo.oligo_name,
            ovo.date_collected
            FROM oligo_variant_overlaps ovo
              JOIN variant_sites ON variant_sites.id = ovo.variant_id
              JOIN fasta_records ON fasta_records.id = variant_sites.fasta_record_id
              JOIN oligos ON oligos.primer_set_id = ovo.primer_set_id AND oligos.name::text = ovo.oligo_name::text
          WHERE (ovo.oligo_end::numeric - ovo.coords) <= 5::numeric
          GROUP BY ovo.detailed_geo_location_id, oligos.locus, ovo.primer_set_id, ovo.oligo_name, ovo.date_collected
        )
        SELECT all_combos.date_collected,
          all_combos.detailed_geo_location_id,
          COALESCE(variants.observed_count, 0::bigint) AS observed_count,
          COALESCE(variants.three_p_score, 0::double precision) AS three_p_score,
          all_combos.locus,
          all_combos.primer_set_id,
          all_combos.oligo_name
          FROM all_combos
            LEFT JOIN variants ON all_combos.date_collected = variants.date_collected AND all_combos.primer_set_id = variants.primer_set_id AND all_combos.oligo_name::text = variants.oligo_name::text AND all_combos.detailed_geo_location_id = variants.detailed_geo_location_id
      WITH DATA;


    CREATE MATERIALIZED VIEW IF NOT EXISTS identify_primers_for_notifications AS
    WITH first_query AS (
                          SELECT primer_set_subscriptions.user_id,
                                 primer_set_subscriptions.primer_set_id,
                                 join_subscribed_location_to_ids.detailed_geo_location_id,
                                 primer_sets.name AS set_name,
                                                     oligos.name AS primer_name,
                                                                    oligos.id as oligo_id,
                                                                                 users.lookback_days,
                                                                                 users.variant_fraction_threshold,
                                                                                 oligo_variant_overlaps.region,
                                                                                 oligo_variant_overlaps.subregion,
                                                                                 oligo_variant_overlaps.division,
                                                                                 oligo_variant_overlaps.subdivision,
                                                                                 oligo_variant_overlaps.detailed_geo_location_id AS unused_id,
                                                                                                                                    oligo_variant_overlaps.coords,
                                                                                                                                    count(oligo_variant_overlaps.variant_id) AS variant_count
    FROM primer_set_subscriptions
    JOIN primer_sets ON primer_sets.id = primer_set_subscriptions.primer_set_id
    JOIN oligos ON primer_sets.id = oligos.primer_set_id
    JOIN oligo_variant_overlaps ON oligo_variant_overlaps.oligo_id = oligos.id
    JOIN join_subscribed_location_to_ids ON join_subscribed_location_to_ids.user_id = primer_set_subscriptions.user_id and join_subscribed_location_to_ids.detailed_geo_location_id = oligo_variant_overlaps.detailed_geo_location_id
    JOIN users ON users.id = primer_set_subscriptions.user_id
    WHERE oligo_variant_overlaps.date_collected >= (CURRENT_DATE - users.lookback_days)
    GROUP BY primer_set_subscriptions.user_id, primer_set_subscriptions.primer_set_id, primer_sets.name, oligos.id, oligos.name,
             join_subscribed_location_to_ids.detailed_geo_location_id,
             users.lookback_days, users.variant_fraction_threshold,
             oligo_variant_overlaps.region, oligo_variant_overlaps.subregion, oligo_variant_overlaps.division, oligo_variant_overlaps.subdivision, oligo_variant_overlaps.coords, oligo_variant_overlaps.detailed_geo_location_id
    ), second_query AS (
                         SELECT fasta_records.detailed_geo_location_id,
                                count(fasta_records.id) AS records_count,
                                                           users.lookback_days
    FROM fasta_records
    JOIN join_subscribed_location_to_ids ON join_subscribed_location_to_ids.detailed_geo_location_id = fasta_records.detailed_geo_location_id
    JOIN users ON join_subscribed_location_to_ids.user_id = users.id
    WHERE fasta_records.date_collected >= (CURRENT_DATE - users.lookback_days)
    GROUP BY fasta_records.detailed_geo_location_id, users.lookback_days
    HAVING count(fasta_records.id) >= 20
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
           first_query.variant_fraction_threshold,
           second_query.detailed_geo_location_id,
           second_query.records_count,
           first_query.variant_count::numeric / second_query.records_count::numeric AS fraction_variant
    FROM first_query
    JOIN second_query ON second_query.detailed_geo_location_id = first_query.detailed_geo_location_id
    WHERE (first_query.variant_count::numeric / second_query.records_count::numeric)::double precision >= first_query.variant_fraction_threshold
    WITH DATA;
    
    grant select on all tables in schema public to primer_monitor_ro;
    SQL
  end
end
