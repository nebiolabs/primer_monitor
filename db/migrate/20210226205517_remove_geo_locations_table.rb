class RemoveGeoLocationsTable < ActiveRecord::Migration[6.1]
  def up
    remove_column :fasta_records, :geo_location_id
    remove_column :subscribed_geo_locations, :geo_location_id
    remove_column :subscribed_geo_locations, :detailed_geo_locations_id
    drop_table :geo_locations
  end

  def down
    create_table :geo_locations do |t|
      t.string :parent_location
      t.string :region, null:false
      t.string :division, null:false

      t.timestamps
    end
    add_reference :fasta_records, :geo_location
    add_reference :subscribed_geo_locations, :geo_location
    add_reference :subscribed_geo_locations, :detailed_geo_locations
  end
end
