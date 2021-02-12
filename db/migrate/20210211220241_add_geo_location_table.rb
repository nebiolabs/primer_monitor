class AddGeoLocationTable < ActiveRecord::Migration[6.1]
  def change
    create_table :geo_location do |t|
      t.string :parent_location
      t.string :region, null:false
      t.string :division, null:false

      t.timestamps
    end
  end
end
