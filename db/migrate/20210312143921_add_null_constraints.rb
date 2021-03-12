class AddNullConstraints < ActiveRecord::Migration[6.1]
  def change
    change_column_null :verified_notifications, :user_id, false
    change_column_null :verified_notifications, :status, false
    change_column_null :primer_set_subscriptions, :user_id, false
    change_column_null :primer_set_subscriptions, :primer_set_id, false
    change_column_null :subscribed_geo_locations, :detailed_geo_location_alias_id, false # Some of these are null still...
    change_column_null :detailed_geo_location_aliases, :world, false
    change_column_null :location_alias_joins, :detailed_geo_location_alias_id, false
    change_column_null :location_alias_joins, :detailed_geo_location_id, false
    change_column_null :fasta_records, :strain, false
    change_column_null :fasta_records, :detailed_geo_location_id, false
    change_column_null :variant_sites, :fasta_record_id, false
    change_column_null :variant_sites, :ref_start, false
    change_column_null :variant_sites, :ref_end, false
    change_column_null :variant_sites, :variant, false
    change_column_null :variant_sites, :variant_type, false
  end
end
