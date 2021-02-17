class RenameGeoLocationTable < ActiveRecord::Migration[6.1]
  def change
    rename_table :geo_location, :geo_locations
  end
end
