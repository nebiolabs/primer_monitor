class AddGeoLocationIdIndex < ActiveRecord::Migration[6.1]
  def change
    add_index :fasta_records, :geo_location_id
    add_foreign_key :fasta_records, :geo_locations
  end
end
