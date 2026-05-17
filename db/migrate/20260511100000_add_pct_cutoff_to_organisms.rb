# frozen_string_literal: true

class AddPctCutoffToOrganisms < ActiveRecord::Migration[7.0]
  def change
    add_column :organisms, :pct_cutoff, :float, null: false, default: 1.0
  end
end
