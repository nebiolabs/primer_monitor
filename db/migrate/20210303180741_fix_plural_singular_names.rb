class FixPluralSingularNames < ActiveRecord::Migration[6.1]
  def up
    rename_column :proposed_notifications, :primer_sets_id, :primer_set_id
    rename_column :proposed_notifications, :oligos_id, :oligo_id
    rename_column :proposed_notifications, :verified_notifications_id, :verified_notification_id
    rename_column :verified_notifications, :users_id, :user_id
    rename_column :location_alias_join, :detailed_geo_locations_id, :detailed_geo_location_id
    rename_column :location_alias_join, :detailed_geo_location_aliases_id, :detailed_geo_location_alias_id

    rename_table :location_alias_join, :location_alias_joins

    execute <<-SQL
    ALTER MATERIALIZED VIEW identify_primers_for_notification RENAME TO identify_primers_for_notifications;
    ALTER VIEW join_subscribed_location_to_id RENAME TO join_subscribed_location_to_ids;
    SQL
  end
  
  def down 
    rename_column :proposed_notifications, :primer_set_id, :primer_sets_id
    rename_column :proposed_notifications, :oligo_id, :oligos_id
    rename_column :proposed_notifications, :verified_notification_id, :verified_notifications_id
    rename_column :verified_notifications, :user_id, :users_id
    rename_column :location_alias_joins, :detailed_geo_location_id, :detailed_geo_locations_id
    rename_column :location_alias_joins, :detailed_geo_location_alias_id, :detailed_geo_location_aliases_id

    rename_table :location_alias_joins, :location_alias_join

    execute <<-SQL
    ALTER MATERIALIZED VIEW identify_primers_for_notifications RENAME TO identify_primers_for_notification;
    ALTER VIEW join_subscribed_location_to_ids RENAME TO join_subscribed_location_to_id;
    SQL
  end
end
