# frozen_string_literal: true

class AddAuthlogicDetailsToUser < ActiveRecord::Migration[6.1]
  def change
    change_table :users do |t|
      t.string :login

      # Authlogic::ActsAsAuthentic::Password
      t.string :crypted_password
      t.string :password_salt

      # Authlogic::ActsAsAuthentic::PersistenceToken
      t.string :persistence_token
      t.index :persistence_token, unique: true

      # Authlogic::ActsAsAuthentic::SingleAccessToken
      t.string :single_access_token
      t.index :single_access_token, unique: true

      # Authlogic::ActsAsAuthentic::PerishableToken
      t.string :perishable_token
      t.index :perishable_token, unique: true

      # Authlogic auto managed fields
      t.integer :login_count, default: 0, null: false
      t.integer :failed_login_count, default: 0, null: false
      t.datetime :last_request_at
      t.datetime :current_login_at
      t.datetime :last_login_at
      t.string :current_login_ip
      t.string :last_login_ip

      # See "Magic States" in Authlogic::Session::Base
      t.boolean :active, default: false
      t.boolean :approved, default: false
      t.boolean :confirmed, default: false

    end
  end
end
