class CreateTableJoiningLocationAndAlias < ActiveRecord::Migration[6.1]
  def change
    create_table :location_alias_join do |t|
      t.references :detailed_geo_locations, foreign_key: true
      t.references :detailed_geo_location_aliases, foreign_key: true, index: false

      t.timestamps
    end
    add_index :location_alias_join, :detailed_geo_locations_id, unique: true, name: 'alias_join_index'
  end
end
