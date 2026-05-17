# frozen_string_literal: true

class CreateLineageVariantPrimerOverlaps < ActiveRecord::Migration[7.0]
  def change
    create_table :lineage_variant_primer_overlaps do |t|
      t.references :organism, null: false, foreign_key: true
      t.string :lineage_group_key, null: false
      t.integer :ref_start, null: false
      t.integer :ref_end, null: false
      t.string :variant_type, null: false
      t.string :variant, null: false
      t.float :frequency_pct, null: false
      t.references :oligo, null: false, foreign_key: true
      t.datetime :created_at, null: false, default: -> { 'NOW()' }
    end

    add_index :lineage_variant_primer_overlaps,
              %i[organism_id lineage_group_key],
              name: 'lvpo_organism_lineage_group_key'

    add_index :lineage_variant_primer_overlaps,
              %i[organism_id lineage_group_key ref_start variant_type variant oligo_id],
              unique: true,
              name: 'lineage_variant_primer_overlaps_unique'
  end
end
