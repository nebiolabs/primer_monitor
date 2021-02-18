# frozen_string_literal: true

class AddDeviseToUsers < ActiveRecord::Migration[6.1]
  def self.up
    add_column :users, :reset_password_token, :string
    add_column :users, :reset_password_sent_at, :datetime
    add_column :users, :remember_token, :string
    add_column :users, :remember_created_at, :datetime
    add_column :users, :authentication_token, :string
    add_column :users, :confirmation_token, :string, limit: 255
    add_column :users, :confirmed_at, :timestamp
    add_column :users, :confirmation_sent_at, :timestamp
    add_column :users, :unconfirmed_email, :string

    ## Lockable
    add_column :users, :failed_attempts, :integer, default: 0, null: false
    add_column :users, :unlock_token, :string
    add_column :users, :locked_at, :datetime

    execute "UPDATE users SET confirmed_at = created_at, confirmation_sent_at = created_at"
    remove_index :users, :perishable_token
    remove_index :users, :persistence_token
    remove_index :users, :single_access_token

    add_index :users, :email,                unique: true
    add_index :users, :unconfirmed_email,    unique: true
    add_index :users, :reset_password_token, unique: true
    add_index :users, :confirmation_token,   unique: true
    add_index :users, :unlock_token,         unique: true


    rename_column :users, :crypted_password, :encrypted_password
    ## Trackable
    rename_column :users, :current_login_at, :current_sign_in_at
    rename_column :users, :last_login_at, :last_sign_in_at
    rename_column :users, :current_login_ip, :current_sign_in_ip
    rename_column :users, :last_login_ip, :last_sign_in_ip
    rename_column :users, :login_count, :sign_in_count

    remove_column :users, :persistence_token
    remove_column :users, :single_access_token
    remove_column :users, :perishable_token
    remove_column :users, :password_salt
  end

  def self.down


    remove_index :users, :email
    remove_index :users, :unconfirmed_email
    remove_index :users, :reset_password_token
    remove_index :users, :confirmation_token
    remove_index :users, :unlock_token


    add_column :users, :perishable_token, :string
    add_column :users, :single_access_token, :string
    add_column :users, :persistence_token, :string
    add_column :users, :password_salt, :string

    add_index :users, :perishable_token
    add_index :users, :persistence_token
    add_index :users, :single_access_token

    rename_column :users, :encrypted_password, :crypted_password
    rename_column :users, :current_sign_in_at, :current_login_at
    rename_column :users, :last_sign_in_at, :last_login_at
    rename_column :users, :current_sign_in_ip, :current_login_ip
    rename_column :users, :last_sign_in_ip, :last_login_ip
    rename_column :users, :sign_in_count, :login_count

    remove_column :users, :confirmation_token
    remove_column :users, :confirmed_at
    remove_column :users, :confirmation_sent_at
    remove_column :users, :authentication_token
    remove_column :users, :remember_created_at
    remove_column :users, :remember_token
    remove_column :users, :reset_password_sent_at
    remove_column :users, :reset_password_token

    remove_column :users, :failed_attempts
    remove_column :users, :unlock_token
    remove_column :users, :locked_at
    remove_column :users, :unconfirmed_email
  end
end
