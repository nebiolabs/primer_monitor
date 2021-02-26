class CreateGeoLocationAliasTable < ActiveRecord::Migration[6.1]
  def change
    create_table :detailed_geo_location_aliases do |t|
      t.string :world
      t.string :region
      t.string :subregion
      t.string :division
      t.string :subdivision
      t.string :locality
      t.string :sublocality
      t.float :latitude
      t.float :longitude

      t.timestamps
    end
    add_index :detailed_geo_location_aliases, [:region, :subregion, :division, :subdivision, :locality, :sublocality, :latitude, :longitude], unique: true, name: 'alias_full_record'

    add_reference :subscribed_geo_locations, :detailed_geo_location_alias, foreign_key: true, index: false
    add_index :subscribed_geo_locations, :detailed_geo_location_alias_id, name: 'tmp'

  end
end
