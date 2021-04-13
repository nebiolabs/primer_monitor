class CreateInitialScoreMatView < ActiveRecord::Migration[6.1]
  def up
    execute "
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
    "
  end

  def down
    execute <<-SQL
    DROP MATERIALIZED VIEW IF EXISTS initial_score;
    SQL
  end
end
