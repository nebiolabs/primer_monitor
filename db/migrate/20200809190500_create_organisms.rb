# frozen_string_literal: true

class CreateOrganisms < ActiveRecord::Migration[6.0]
  def change
    create_table :organisms do |t|
      t.integer :ncbi_taxon_id, null: false, unique: true
      t.string :name, null: false, unique: true

      t.timestamps
    end
  end
end
