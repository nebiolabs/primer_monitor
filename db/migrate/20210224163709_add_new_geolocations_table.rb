class AddNewGeolocationsTable < ActiveRecord::Migration[6.1]
  def change
    create_table :detailed_geo_locations do |t|
      t.string :world, null: false
      t.string :region, null: false
      t.string :subregion
      t.string :division
      t.string :subdivision
      t.string :locality
      t.float :latitude, null: false
      t.float :longitude, null: false

      t.timestamps
    end
    add_index :detailed_geo_locations, [:region, :subregion, :division, :subdivision, :locality, :latitude, :longitude], unique: true, name: 'full_record'

    add_reference :fasta_records, :detailed_geo_locations, foreign_key: true
    add_reference :subscribed_geo_locations, :detailed_geo_locations, foreign_key: true

  end
end
