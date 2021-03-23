class JoinLocationSubscriptionToId < ActiveRecord::Migration[6.1]
  def up
    execute "
    CREATE VIEW join_subscribed_location_to_id AS
      WITH subscribed_ids AS (
        SELECT subscribed_geo_locations.user_id, region, subregion, division, subdivision FROM detailed_geo_location_aliases
          INNER JOIN subscribed_geo_locations ON subscribed_geo_locations.detailed_geo_location_alias_id = detailed_geo_location_aliases.id
      )
      SELECT subscribed_ids.user_id, detailed_geo_locations.id AS detailed_geo_location_id FROM detailed_geo_locations
        INNER JOIN subscribed_ids ON (subscribed_ids.region IS NULL OR subscribed_ids.region = detailed_geo_locations.region) AND 
            (subscribed_ids.subregion IS NULL OR subscribed_ids.subregion = detailed_geo_locations.subregion) AND 
            (subscribed_ids.division IS NULL OR subscribed_ids.division = detailed_geo_locations.division) AND
            (subscribed_ids.subdivision IS NULL OR subscribed_ids.subdivision = detailed_geo_locations.subdivision)
    "
  end

  def down
    execute <<-SQL
      DROP VIEW IF EXISTS join_subscribed_location_to_id;
    SQL
  end
end
