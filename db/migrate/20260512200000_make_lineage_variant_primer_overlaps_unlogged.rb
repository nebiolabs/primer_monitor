# frozen_string_literal: true

class MakeLineageVariantPrimerOverlapsUnlogged < ActiveRecord::Migration[7.1]
  def up
    execute 'ALTER TABLE lineage_variant_primer_overlaps SET UNLOGGED'
  end

  def down
    execute 'ALTER TABLE lineage_variant_primer_overlaps SET LOGGED'
  end
end
