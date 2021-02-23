class CreateSubscribedGeoLocations < ActiveRecord::Migration[6.1]

  def up
    create_table :subscribed_geo_locations do |t|
      t.references :user, foreign_key: true, null: false
      t.references :geo_location, foreign_key: true, null: false
      t.timestamps
    end
    add_index :subscribed_geo_locations, [:user_id, :geo_location_id], unique: true
    world = GeoLocation.find_or_create_by(region: 'world') do |loc|
      loc.division = 'world'
    end
    User.all.each do |user|
      SubscribedGeoLocation.find_or_create_by(user_id: user.id, geo_location_id: world.id )
    end
  end

  def down
    remove_index :subscribed_geo_locations, [:user_id, :geo_location_id]
    drop_table :subscribed_geo_locations

    GeoLocation.update_all(parent_location: nil)
    GeoLocation.find_by(region: 'world')&.destroy
  end
end
