class RenameDetailedGeoLocationsIdColumn < ActiveRecord::Migration[6.1]
  def change
    rename_column :fasta_records, :detailed_geo_locations_id, :detailed_geo_location_id
  end
end
