class RemoveLocationAliasJoin < ActiveRecord::Migration[6.1]
  def up
    add_reference :detailed_geo_locations, :detailed_geo_location_alias, foreign_key: true
    execute <<~SQL.chop
      update detailed_geo_locations 
      set detailed_geo_location_alias_id = location_alias_joins.detailed_geo_location_alias_id
      from location_alias_joins
      where location_alias_joins.detailed_geo_location_id = detailed_geo_locations.id
    SQL

    change_column_null :detailed_geo_locations, :detailed_geo_location_alias_id, false
    drop_table :location_alias_joins
  end

  def down
    create_table :location_alias_joins do |t|
      t.references :detailed_geo_location, foreign_key: true, null:false
      t.references :detailed_geo_location_alias, foreign_key: true, index: false, null:false

      t.timestamps
    end

    execute <<~SQL.chop
      insert into location_alias_joins (detailed_geo_location_id, detailed_geo_location_alias_id, created_at, updated_at)
      select  id, detailed_geo_location_alias_id, created_at, updated_at  
      from detailed_geo_locations
    SQL
    add_index :location_alias_joins, :detailed_geo_location_id, unique: true, name: 'alias_join_index'

    remove_column :detailed_geo_locations, :detailed_geo_location_alias_id

  end
end
