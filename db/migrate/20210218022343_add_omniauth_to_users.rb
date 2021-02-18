class AddOmniauthToUsers < ActiveRecord::Migration[6.1]
  def change
    add_column :users, :provider, :string
    add_column :users, :uid, :string

    remove_index :users, :openid_identifier
    remove_column :users, :openid_identifier
  end
end
