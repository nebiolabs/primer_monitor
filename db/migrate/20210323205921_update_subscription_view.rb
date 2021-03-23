class UpdateSubscriptionView < ActiveRecord::Migration[6.1]
  def up
    execute "
    DROP VIEW IF EXISTS join_subscribed_location_to_id;

    CREATE VIEW join_subscribed_location_to_id AS
    WITH subscribed_ids AS (
      SELECT subscribed_geo_locations.user_id,
        subscribed_geo_locations.detailed_geo_location_alias_id,
        detailed_geo_location_aliases.region,
        detailed_geo_location_aliases.subregion,
        detailed_geo_location_aliases.division,
        detailed_geo_location_aliases.subdivision
        FROM detailed_geo_location_aliases
          JOIN subscribed_geo_locations ON subscribed_geo_locations.detailed_geo_location_alias_id = detailed_geo_location_aliases.id
    )
    SELECT subscribed_ids.user_id, subscribed_ids.detailed_geo_location_alias_id as subscribed_id,
    detailed_geo_location_aliases.id AS detailed_geo_location_alias_id,
    detailed_geo_locations.id as detailed_geo_location_id
    FROM detailed_geo_location_aliases
    JOIN subscribed_ids ON (subscribed_ids.region IS NULL OR subscribed_ids.region::text = detailed_geo_location_aliases.region::text) AND (subscribed_ids.subregion IS NULL OR subscribed_ids.subregion::text = detailed_geo_location_aliases.subregion::text) AND (subscribed_ids.division IS NULL OR subscribed_ids.division::text = detailed_geo_location_aliases.division::text) AND (subscribed_ids.subdivision IS NULL OR subscribed_ids.subdivision::text = detailed_geo_location_aliases.subdivision::text)     
    join detailed_geo_locations on detailed_geo_locations.detailed_geo_location_alias_id = detailed_geo_location_aliases.id;
    "
  end

  def down
    execute <<-SQL
      DROP VIEW IF EXISTS join_subscribed_location_to_id;
    SQL
  end
end