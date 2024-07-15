class AddMissingUniqueIndices < ActiveRecord::Migration[6.1]
  def change
    add_index :roles, :name, unique: true
    # abbreviated name due to index name too long error
    add_index :subscribed_geo_locations, %i[user_id detailed_geo_location_alias_id], unique: true, name: 'index_subscribed_geo_locs_on_user_and_detailed_geo_loc_alias_id'
    add_index :user_roles, %i[user_id role_id], unique: true
  end
end
