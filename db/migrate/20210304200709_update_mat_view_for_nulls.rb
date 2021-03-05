class UpdateMatViewForNulls < ActiveRecord::Migration[6.1]
  def up
    execute "
      DROP MATERIALIZED VIEW IF EXISTS identify_primers_for_notifications;
      DROP MATERIALIZED VIEW IF EXISTS oligo_variant_overlaps;    

      CREATE MATERIALIZED VIEW IF NOT EXISTS variant_overlaps AS
        SELECT oligos.id AS oligo_id,
          oligos.name AS oligo_name,
          oligos.ref_start AS oligo_start,
          oligos.ref_end AS oligo_end,
          oligos.short_name as oligo_short_name,
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
        JOIN oligos ON variant_sites.ref_start >= oligos.ref_start AND variant_sites.ref_start <= oligos.ref_end OR variant_sites.ref_end >= oligos.ref_start AND variant_sites.ref_end <= oligos.ref_end OR variant_sites.ref_start < oligos.ref_start AND variant_sites.ref_end > oligos.ref_end
        JOIN primer_sets on oligos.primer_set_id = primer_sets.id
        WHERE (variant_sites.variant_type::text = 'D'::text OR variant_sites.variant_type::text = 'X'::text OR variant_sites.variant_type::text = 'I'::text) AND variant_sites.variant::text !~~ '%N%'::text
      WITH DATA;

      CREATE INDEX ON variant_overlaps(region, subregion, division, subdivision, date_collected);

      CREATE MATERIALIZED VIEW IF NOT EXISTS time_counts AS
        WITH region_time_count AS (
          SELECT COUNT(*) AS region_time_count, COALESCE(detailed_geo_locations.region, '') as region, coalesce(fasta_records.date_collected, '1900-01-01') as date_collected FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, fasta_records.date_collected     
          ),
          region_subregion_time_count AS (
          SELECT COUNT(*) AS region_subregion_time_count, COALESCE(detailed_geo_locations.region, '') as region, COALESCE(detailed_geo_locations.subregion, '') as subregion, coalesce(fasta_records.date_collected, '1900-01-01') as date_collected FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, fasta_records.date_collected
          ),
          region_subregion_division_time_count AS (
          SELECT COUNT(*) AS region_subregion_division_time_count, COALESCE(detailed_geo_locations.region, '') as region, COALESCE(detailed_geo_locations.subregion, '') as subregion, COALESCE(detailed_geo_locations.division, '') as division, coalesce(fasta_records.date_collected, '1900-01-01') as date_collected FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, fasta_records.date_collected
          ),
          region_subregion_division_subdivision_time_count AS (
          SELECT COUNT(*) AS region_subregion_division_subdivision_time_count, COALESCE(detailed_geo_locations.region, '') as region, COALESCE(detailed_geo_locations.subregion, '') as subregion, COALESCE(detailed_geo_locations.division, '') as division, COALESCE(detailed_geo_locations.subdivision, '') as subdivision, coalesce(fasta_records.date_collected, '1900-01-01') as date_collected FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision, fasta_records.date_collected
          )
        SELECT 
          region_subregion_division_subdivision_time_count.region, region_subregion_division_subdivision_time_count.subregion, region_subregion_division_subdivision_time_count.division, region_subregion_division_subdivision_time_count.subdivision, region_subregion_division_subdivision_time_count.date_collected,
          region_time_count.region_time_count, region_subregion_time_count.region_subregion_time_count, region_subregion_division_time_count.region_subregion_division_time_count,  region_subregion_division_subdivision_time_count.region_subregion_division_subdivision_time_count
        FROM region_subregion_division_subdivision_time_count
        INNER JOIN region_time_count ON region_time_count.region = region_subregion_division_subdivision_time_count.region and region_time_count.date_collected = region_subregion_division_subdivision_time_count.date_collected
        INNER JOIN region_subregion_time_count ON region_subregion_time_count.region = region_subregion_division_subdivision_time_count.region AND region_subregion_time_count.subregion = region_subregion_division_subdivision_time_count.subregion and region_subregion_time_count.date_collected = region_subregion_division_subdivision_time_count.date_collected
        INNER JOIN region_subregion_division_time_count on region_subregion_division_time_count.region = region_subregion_division_subdivision_time_count.region AND region_subregion_division_time_count.subregion = region_subregion_division_subdivision_time_count.subregion and region_subregion_division_time_count.division = region_subregion_division_subdivision_time_count.division and region_subregion_division_time_count.date_collected = region_subregion_division_subdivision_time_count.date_collected
      WITH DATA;

      CREATE INDEX ON time_counts(region, subregion, division, subdivision, date_collected);

      CREATE MATERIALIZED VIEW IF NOT EXISTS counts AS
        WITH region_count AS (
          SELECT COUNT(*) AS region_count, COALESCE(detailed_geo_locations.region, '') as region FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region       
          ),
          region_subregion_count AS (
          SELECT COUNT(*) AS region_subregion_count, COALESCE(detailed_geo_locations.region, '') as region, COALESCE(detailed_geo_locations.subregion, '') as subregion FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion
          ),
          region_subregion_division_count AS (
          SELECT COUNT(*) AS region_subregion_division_count, COALESCE(detailed_geo_locations.region, '') as region, COALESCE(detailed_geo_locations.subregion, '') as subregion, COALESCE(detailed_geo_locations.division, '') as division FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division
          ),
          region_subregion_division_subdivision_count AS (
          SELECT COUNT(*) AS region_subregion_division_subdivision_count, COALESCE(detailed_geo_locations.region, '') as region, COALESCE(detailed_geo_locations.subregion, '') as subregion, COALESCE(detailed_geo_locations.division, '') as division, COALESCE(detailed_geo_locations.subdivision, '') as subdivision FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision
          )
        SELECT 
          region_subregion_division_subdivision_count.region, region_subregion_division_subdivision_count.subregion, region_subregion_division_subdivision_count.division, region_subregion_division_subdivision_count.subdivision,
          region_count.region_count, region_subregion_count.region_subregion_count, region_subregion_division_count.region_subregion_division_count, region_subregion_division_subdivision_count.region_subregion_division_subdivision_count
        FROM region_subregion_division_subdivision_count
        INNER JOIN region_count on region_subregion_division_subdivision_count.region = region_count.region
        INNER JOIN region_subregion_count ON region_subregion_count.region = region_subregion_division_subdivision_count.region AND region_subregion_count.subregion = region_subregion_division_subdivision_count.subregion
        INNER JOIN region_subregion_division_count on region_subregion_division_count.region = region_subregion_division_subdivision_count.region AND region_subregion_division_count.subregion = region_subregion_division_subdivision_count.subregion AND region_subregion_division_count.division = region_subregion_division_subdivision_count.division
      WITH DATA;

      CREATE INDEX ON counts(region, subregion, division, subdivision);

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

      CREATE MATERIALIZED VIEW IF NOT EXISTS identify_primers_for_notifications AS
        WITH first_query AS (
         SELECT primer_set_subscriptions.user_id,
            primer_set_subscriptions.primer_set_id,
            join_subscribed_location_to_ids.detailed_geo_location_id,
            primer_sets.name AS set_name,
            oligos.name AS primer_name,
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
          GROUP BY primer_set_subscriptions.user_id, primer_set_subscriptions.primer_set_id, primer_sets.name, oligos.name,
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
      "
  end

  def down
    execute <<-SQL
    DROP MATERIALIZED VIEW IF EXISTS identify_primers_for_notifications;
    DROP MATERIALIZED VIEW IF EXISTS oligo_variant_overlaps;
    DROP MATERIALIZED VIEW IF EXISTS variant_overlaps;
    DROP MATERIALIZED VIEW IF EXISTS time_counts;
    DROP MATERIALIZED VIEW IF EXISTS counts;
    SQL
  end
end