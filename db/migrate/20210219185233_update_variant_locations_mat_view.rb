class UpdateVariantLocationsMatView < ActiveRecord::Migration[6.1]
  def up
    execute "
      DROP MATERIALIZED VIEW IF EXISTS oligo_variant_overlaps;

      CREATE MATERIALIZED VIEW IF NOT EXISTS oligo_variant_overlaps AS
      WITH big_query AS (
        SELECT oligos.id AS oligo_id, oligos.name AS oligo_name, oligos.ref_start AS oligo_start, oligos.ref_end AS oligo_end,
              variant_sites.id AS variant_id, variant_sites.variant_type, variant_sites.variant, variant_sites.ref_start AS variant_start, variant_sites.ref_end AS variant_end,
              geo_locations.region, geo_locations.division, geo_locations.id as geo_id,
              fasta_records.date_collected FROM variant_sites
        INNER JOIN fasta_records ON variant_sites.fasta_record_id = fasta_records.id
        INNER JOIN geo_locations ON fasta_records.geo_location_id = geo_locations.id
        INNER JOIN oligos ON (variant_sites.ref_start BETWEEN oligos.ref_start AND oligos.ref_end)
                  OR (variant_sites.ref_end BETWEEN oligos.ref_start AND oligos.ref_end)
                  OR (variant_sites.ref_start < oligos.ref_start AND variant_sites.ref_end > oligos.ref_end)
        WHERE (variant_type = 'D' OR variant_type = 'X') AND variant NOT LIKE '%N%'
        ),
        insert_query as (
        SELECT oligos.id AS oligo_id, oligos.name AS oligo_name, oligos.ref_start AS oligo_start, oligos.ref_end AS oligo_end,
              variant_sites.id AS variant_id, variant_sites.variant_type, variant_sites.variant, variant_sites.ref_start AS variant_start, variant_sites.ref_end AS variant_end,
              geo_locations.region, geo_locations.division, geo_locations.id as geo_id,
              fasta_records.date_collected FROM variant_sites
        INNER JOIN fasta_records ON variant_sites.fasta_record_id = fasta_records.id
        INNER JOIN geo_locations ON fasta_records.geo_location_id = geo_locations.id
        INNER JOIN oligos ON (variant_sites.ref_start BETWEEN oligos.ref_start AND oligos.ref_end)
                  OR (variant_sites.ref_end BETWEEN oligos.ref_start AND oligos.ref_end)
                  OR (variant_sites.ref_start < oligos.ref_start AND variant_sites.ref_end > oligos.ref_end)
        WHERE variant_type = 'I' AND variant NOT LIKE '%N%'
        ),
        region_count AS (
        SELECT COUNT(*) AS region_count, geo_locations.region FROM fasta_records 
          INNER JOIN geo_locations ON fasta_records.geo_location_id = geo_locations.id
          GROUP BY geo_locations.region       
        ),
        region_division_count AS (
        SELECT COUNT(*) AS region_division_count, geo_locations.region, geo_locations.division FROM fasta_records 
          INNER JOIN geo_locations ON fasta_records.geo_location_id = geo_locations.id
          GROUP BY geo_locations.id, geo_locations.region, geo_locations.division
        ),
        region_time_count AS (
        SELECT COUNT(*) AS region_time_count, geo_locations.region, fasta_records.date_collected FROM fasta_records 
          INNER JOIN geo_locations ON fasta_records.geo_location_id = geo_locations.id
          GROUP BY geo_locations.region, fasta_records.date_collected    
        ),
        region_division_time_count AS (
        SELECT COUNT(*) AS region_division_time_count, geo_locations.region, geo_locations.division, fasta_records.date_collected FROM fasta_records 
          INNER JOIN geo_locations ON fasta_records.geo_location_id = geo_locations.id
          GROUP BY geo_locations.region, geo_locations.division, fasta_records.date_collected    
        )
        SELECT big_query.*, region_count.region_count, region_division_count.region_division_count, region_time_count.region_time_count, region_division_time_count.region_division_time_count,
          GENERATE_SERIES(
            LOWER(numrange(coord_overlaps.oligo_start, coord_overlaps.oligo_end) * numrange(coord_overlaps.variant_start, coord_overlaps.variant_end)),
            UPPER(numrange(coord_overlaps.oligo_start, coord_overlaps.oligo_end) * numrange(coord_overlaps.variant_start, coord_overlaps.variant_end)) - 1)
          AS coords FROM big_query AS coord_overlaps
        INNER JOIN big_query ON coord_overlaps.oligo_id = big_query.oligo_id AND coord_overlaps.variant_id = big_query.variant_id
        INNER JOIN region_count ON region_count.region = big_query.region
        INNER JOIN region_division_count on region_division_count.region = big_query.region AND region_division_count.division = big_query.division
        INNER JOIN region_time_count ON region_time_count.region = big_query.region AND region_time_count.date_collected = big_query.date_collected
        INNER JOIN region_division_time_count ON region_division_time_count.region = big_query.region AND region_division_time_count.division = big_query.division AND region_division_time_count.date_collected = big_query.date_collected
      UNION
        SELECT insert_query.*, region_count.region_count, region_division_count.region_division_count, region_time_count.region_time_count, region_division_time_count.region_division_time_count,
          insert_query.variant_start AS coords FROM insert_query AS coord_overlaps
        INNER JOIN insert_query ON coord_overlaps.oligo_id = insert_query.oligo_id AND coord_overlaps.variant_id = insert_query.variant_id
        INNER JOIN region_count ON region_count.region = insert_query.region
        INNER JOIN region_division_count ON region_division_count.region = insert_query.region AND region_division_count.division = insert_query.division
        INNER JOIN region_time_count ON region_time_count.region = insert_query.region AND region_time_count.date_collected = insert_query.date_collected
        INNER JOIN region_division_time_count ON region_division_time_count.region = insert_query.region AND region_division_time_count.division = insert_query.division AND region_division_time_count.date_collected = insert_query.date_collected
      WITH DATA;
      "
  end

  def down
    execute <<-SQL
      DROP MATERIALIZED VIEW IF EXISTS oligo_variant_overlaps;
    SQL
  end
end
