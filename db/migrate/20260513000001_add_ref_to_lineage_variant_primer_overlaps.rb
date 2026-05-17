class AddRefToLineageVariantPrimerOverlaps < ActiveRecord::Migration[7.1]
  def change
    add_column :lineage_variant_primer_overlaps, :ref, :string
  end
end
