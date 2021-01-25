# frozen_string_literal: true

class CreateUserRoles < ActiveRecord::Migration[6.1]
  def change
    create_table :user_roles do |t|
      t.references :user, foreign_key: true, null: false, index: true
      t.references :role, foreign_key: true, null: false, index: true

      t.timestamps
    end
  end
end
