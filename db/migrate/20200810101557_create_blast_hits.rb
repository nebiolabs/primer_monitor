# frozen_string_literal: true

class CreateBlastHits < ActiveRecord::Migration[6.0]
  def change
    create_table :blast_hits do |t|
      t.references :oligo, foreign_key: true
      t.references :organism, foreign_key: true
      t.integer :num_identities

      t.timestamps
    end
  end
end
