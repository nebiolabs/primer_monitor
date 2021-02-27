class IdentifyPrimersForNotification < ActiveRecord::Migration[6.1]
  def up
    execute "
      DROP MATERIALIZED VIEW IF EXISTS identify_primers_for_notification;

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
    SQL
  end
end
