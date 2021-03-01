class ProperMigrationFormat < ActiveRecord::Migration[6.1]
  def up
    execute "
      DROP MATERIALIZED VIEW IF EXISTS identify_primers_for_notification;
      DROP MATERIALIZED VIEW IF EXISTS oligo_variant_overlaps;

      CREATE MATERIALIZED VIEW IF NOT EXISTS oligo_variant_overlaps AS
      WITH big_query AS (
        SELECT oligos.id AS oligo_id, oligos.name AS oligo_name, oligos.ref_start AS oligo_start, oligos.ref_end AS oligo_end,
              variant_sites.id AS variant_id, variant_sites.variant_type, variant_sites.variant, variant_sites.ref_start AS variant_start, variant_sites.ref_end AS variant_end,
              detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision, detailed_geo_locations.id as detailed_geo_location_id,
              fasta_records.date_collected FROM variant_sites
        INNER JOIN fasta_records ON variant_sites.fasta_record_id = fasta_records.id
        INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
        INNER JOIN oligos ON (variant_sites.ref_start BETWEEN oligos.ref_start AND oligos.ref_end)
                  OR (variant_sites.ref_end BETWEEN oligos.ref_start AND oligos.ref_end)
                  OR (variant_sites.ref_start < oligos.ref_start AND variant_sites.ref_end > oligos.ref_end)
        WHERE (variant_type = 'D' OR variant_type = 'X') AND variant NOT LIKE '%N%'
        ),
        insert_query as (
        SELECT oligos.id AS oligo_id, oligos.name AS oligo_name, oligos.ref_start AS oligo_start, oligos.ref_end AS oligo_end,
              variant_sites.id AS variant_id, variant_sites.variant_type, variant_sites.variant, variant_sites.ref_start AS variant_start, variant_sites.ref_end AS variant_end,
              detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision, detailed_geo_locations.id as detailed_geo_location_id,
              fasta_records.date_collected FROM variant_sites
        INNER JOIN fasta_records ON variant_sites.fasta_record_id = fasta_records.id
        INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
        INNER JOIN oligos ON (variant_sites.ref_start BETWEEN oligos.ref_start AND oligos.ref_end)
                  OR (variant_sites.ref_end BETWEEN oligos.ref_start AND oligos.ref_end)
                  OR (variant_sites.ref_start < oligos.ref_start AND variant_sites.ref_end > oligos.ref_end)
        WHERE variant_type = 'I' AND variant NOT LIKE '%N%'
        ),
        region_count AS (
        SELECT COUNT(*) AS region_count, detailed_geo_locations.region FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region       
        ),
        region_subregion_count AS (
        SELECT COUNT(*) AS region_subregion_count, detailed_geo_locations.region, detailed_geo_locations.subregion FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion
        ),
        region_subregion_division_count AS (
        SELECT COUNT(*) AS region_subregion_division_count, detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division
        ),
        region_subregion_division_subdivision_count AS (
        SELECT COUNT(*) AS region_subregion_division_subdivision_count, detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision
        ),
        region_time_count AS (
        SELECT COUNT(*) AS region_time_count, detailed_geo_locations.region, fasta_records.date_collected FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, fasta_records.date_collected     
        ),
        region_subregion_time_count AS (
        SELECT COUNT(*) AS region_subregion_time_count, detailed_geo_locations.region, detailed_geo_locations.subregion, fasta_records.date_collected FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, fasta_records.date_collected
        ),
        region_subregion_division_time_count AS (
        SELECT COUNT(*) AS region_subregion_division_time_count, detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, fasta_records.date_collected FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, fasta_records.date_collected
        ),
        region_subregion_division_subdivision_time_count AS (
        SELECT COUNT(*) AS region_subregion_division_subdivision_time_count, detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision, fasta_records.date_collected FROM fasta_records 
          INNER JOIN detailed_geo_locations ON fasta_records.detailed_geo_location_id = detailed_geo_locations.id
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision, fasta_records.date_collected
        )
        SELECT big_query.*, region_count.region_count, region_subregion_count.region_subregion_count, region_subregion_division_count.region_subregion_division_count, region_subregion_division_subdivision_count.region_subregion_division_subdivision_count,
                            region_time_count.region_time_count, region_subregion_time_count.region_subregion_time_count, region_subregion_division_time_count.region_subregion_division_time_count, region_subregion_division_subdivision_time_count.region_subregion_division_subdivision_time_count,
          GENERATE_SERIES(
            LOWER(numrange(coord_overlaps.oligo_start, coord_overlaps.oligo_end) * numrange(coord_overlaps.variant_start, coord_overlaps.variant_end)),
            UPPER(numrange(coord_overlaps.oligo_start, coord_overlaps.oligo_end) * numrange(coord_overlaps.variant_start, coord_overlaps.variant_end)) - 1)
          AS coords FROM big_query AS coord_overlaps
        INNER JOIN big_query ON coord_overlaps.oligo_id = big_query.oligo_id AND coord_overlaps.variant_id = big_query.variant_id
        INNER JOIN region_count ON region_count.region = big_query.region
        INNER JOIN region_subregion_count ON region_subregion_count.region = big_query.region AND region_subregion_count.subregion = big_query.subregion
        INNER JOIN region_subregion_division_count on region_subregion_division_count.region = big_query.region AND region_subregion_division_count.subregion = big_query.subregion AND region_subregion_division_count.division = big_query.division
        INNER JOIN region_subregion_division_subdivision_count on region_subregion_division_subdivision_count.region = big_query.region AND region_subregion_division_subdivision_count.subregion = big_query.subregion AND region_subregion_division_subdivision_count.division = big_query.division AND region_subregion_division_subdivision_count.subdivision = big_query.subdivision
        INNER JOIN region_time_count ON region_time_count.region = big_query.region and region_time_count.date_collected = big_query.date_collected
        INNER JOIN region_subregion_time_count ON region_subregion_time_count.region = big_query.region AND region_subregion_time_count.subregion = big_query.subregion and region_subregion_time_count.date_collected = big_query.date_collected
        INNER JOIN region_subregion_division_time_count on region_subregion_division_time_count.region = big_query.region AND region_subregion_division_time_count.subregion = big_query.subregion and region_subregion_division_time_count.division = big_query.division and region_subregion_division_time_count.date_collected = big_query.date_collected
        INNER JOIN region_subregion_division_subdivision_time_count on region_subregion_division_subdivision_time_count.region = big_query.region AND region_subregion_division_subdivision_time_count.subregion = big_query.subregion and region_subregion_division_subdivision_time_count.division = big_query.division and region_subregion_division_subdivision_time_count.subdivision = big_query.subdivision AND region_subregion_division_subdivision_time_count.date_collected = big_query.date_collected
      UNION ALL
        SELECT insert_query.*, region_count.region_count, region_subregion_count.region_subregion_count, region_subregion_division_count.region_subregion_division_count, region_subregion_division_subdivision_count.region_subregion_division_subdivision_count,
                               region_time_count.region_time_count, region_subregion_time_count.region_subregion_time_count, region_subregion_division_time_count.region_subregion_division_time_count, region_subregion_division_subdivision_time_count.region_subregion_division_subdivision_time_count,
          insert_query.variant_start AS coords FROM insert_query AS coord_overlaps
        INNER JOIN insert_query ON coord_overlaps.oligo_id = insert_query.oligo_id AND coord_overlaps.variant_id = insert_query.variant_id
        INNER JOIN region_count ON region_count.region = insert_query.region
        INNER JOIN region_subregion_count ON region_subregion_count.region = insert_query.region AND region_subregion_count.subregion = insert_query.subregion
        INNER JOIN region_subregion_division_count on region_subregion_division_count.region = insert_query.region AND region_subregion_division_count.subregion = insert_query.subregion and region_subregion_division_count.division = insert_query.division
        INNER JOIN region_subregion_division_subdivision_count on region_subregion_division_subdivision_count.region = insert_query.region AND region_subregion_division_subdivision_count.subregion = insert_query.subregion and region_subregion_division_subdivision_count.division = insert_query.division and region_subregion_division_subdivision_count.subdivision = insert_query.subdivision
        INNER JOIN region_time_count ON region_time_count.region = insert_query.region and region_time_count.date_collected = insert_query.date_collected
        INNER JOIN region_subregion_time_count ON region_subregion_time_count.region = insert_query.region AND region_subregion_time_count.subregion = insert_query.subregion and region_subregion_time_count.date_collected = insert_query.date_collected
        INNER JOIN region_subregion_division_time_count on region_subregion_division_time_count.region = insert_query.region AND region_subregion_division_time_count.subregion = insert_query.subregion and region_subregion_division_time_count.division = insert_query.division and region_subregion_division_time_count.date_collected = insert_query.date_collected
        INNER JOIN region_subregion_division_subdivision_time_count on region_subregion_division_subdivision_time_count.region = insert_query.region AND region_subregion_division_subdivision_time_count.subregion = insert_query.subregion and region_subregion_division_subdivision_time_count.division = insert_query.division and region_subregion_division_subdivision_time_count.subdivision = insert_query.subdivision and region_subregion_division_subdivision_time_count.date_collected = insert_query.date_collected
      WITH DATA;

      CREATE MATERIALIZED VIEW IF NOT EXISTS identify_primers_for_notification AS
      with first_query as (
        select primer_set_subscriptions.user_id, primer_set_subscriptions.primer_set_id,
             join_subscribed_location_to_id.detailed_geo_location_id as detailed_geo_location_id,
             primer_sets.name as set_name,
             oligos.name as primer_name,
             users.lookback_days, users.variant_fraction_threshold,
             oligo_variant_overlaps.region,
             oligo_variant_overlaps.division,
             oligo_variant_overlaps.detailed_geo_location_id as unused_id,
             oligo_variant_overlaps.coords,
             count(oligo_variant_overlaps.variant_id) as variant_count from primer_set_subscriptions
        inner join primer_sets on primer_sets.id = primer_set_subscriptions.primer_set_id
        inner join join_subscribed_location_to_id on join_subscribed_location_to_id.user_id = primer_set_subscriptions.user_id
        inner join oligos on primer_sets.id = oligos.primer_set_id
        inner join oligo_variant_overlaps on oligo_variant_overlaps.oligo_id = oligos.id
        inner join users on users.id = primer_set_subscriptions.user_id
        where oligo_variant_overlaps.date_collected >= current_date - users.lookback_days 
        group by primer_set_subscriptions.user_id, primer_set_subscriptions.primer_set_id, primer_sets.name, oligos.name,
             join_subscribed_location_to_id.detailed_geo_location_id,
             users.lookback_days, users.variant_fraction_threshold,
             oligo_variant_overlaps.region,
             oligo_variant_overlaps.division,
             oligo_variant_overlaps.coords,
             oligo_variant_overlaps.detailed_geo_location_id 
        ),
        second_query as (
          select fasta_records.detailed_geo_location_id, count(fasta_records.id) as records_count,
               users.lookback_days 
          from fasta_records
          inner join join_subscribed_location_to_id on join_subscribed_location_to_id.detailed_geo_location_id = fasta_records.detailed_geo_location_id
          inner join users on join_subscribed_location_to_id.user_id = users.id
          where fasta_records.date_collected >= current_date - users.lookback_days
          group by fasta_records.detailed_geo_location_id, users.lookback_days 
          having count(fasta_records.id) >= 20
        )
        select first_query.user_id, first_query.primer_set_id, first_query.set_name, first_query.primer_name, first_query.region, first_query.division, first_query.coords, first_query.variant_count, first_query.variant_fraction_threshold,
             second_query.detailed_geo_location_id, second_query.records_count, first_query.variant_count / second_query.records_count::decimal as fraction_variant from first_query
          inner join second_query on second_query.detailed_geo_location_id = first_query.detailed_geo_location_id
        where (first_query.variant_count / second_query.records_count::decimal) >= first_query.variant_fraction_threshold
      WITH DATA;
      "
  end

  def down
    execute <<-SQL
    DROP MATERIALIZED VIEW IF EXISTS identify_primers_for_notification;
    DROP MATERIALIZED VIEW IF EXISTS oligo_variant_overlaps;
    SQL
  end
end
