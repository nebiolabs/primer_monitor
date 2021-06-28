class AddLocationId < ActiveRecord::Migration[6.1]
  def change
    add_reference :proposed_notifications, :detailed_geo_location_alias, foreign_key: true
  end
end
