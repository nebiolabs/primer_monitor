class RenameGeoLocColumn < ActiveRecord::Migration[6.1]
  def change
    rename_column :fasta_records, :geo_loc_id, :geo_location_id
  end
end
