# frozen_string_literal: true

class AddSeenColumnsToLineageVariantPrimerOverlaps < ActiveRecord::Migration[7.0]
  def change
    add_column :lineage_variant_primer_overlaps, :first_seen_date,     :date
    add_column :lineage_variant_primer_overlaps, :first_seen_lineage,  :string
    add_column :lineage_variant_primer_overlaps, :first_seen_location, :string
    add_column :lineage_variant_primer_overlaps, :last_seen_date,      :date
    add_column :lineage_variant_primer_overlaps, :last_seen_lineage,   :string
    add_column :lineage_variant_primer_overlaps, :last_seen_location,  :string
  end
end
