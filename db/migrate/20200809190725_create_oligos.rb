# frozen_string_literal: true

class CreateOligos < ActiveRecord::Migration[6.0]
  def change
    create_table :oligos do |t|
      t.string :name, null:false
      t.string :sequence, null:false
      t.references :amplicon, foreign_key: true, null:false
      t.bigint :ref_start
      t.bigint :ref_end

      t.timestamps
    end
  end
end
