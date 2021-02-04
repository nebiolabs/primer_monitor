# frozen_string_literal: true

class CreateAmplicons < ActiveRecord::Migration[6.0]
  def change
    create_table :amplicons do |t|
      t.string :name
      t.references :user, foreign_key: true
      t.references :organism, foreign_key: true
      t.timestamps
    end
  end
end
