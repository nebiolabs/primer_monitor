class AddGeoLocIdColumn < ActiveRecord::Migration[6.1]
  def change
    add_column :fasta_records, :geo_loc_id, :integer
  end
end
