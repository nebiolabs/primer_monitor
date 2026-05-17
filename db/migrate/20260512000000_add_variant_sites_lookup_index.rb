# frozen_string_literal: true

class AddVariantSitesLookupIndex < ActiveRecord::Migration[7.0]
  disable_ddl_transaction!

  def change
    add_index :variant_sites, %i[ref_start variant_type variant organism_taxon_id],
              name: 'idx_variant_sites_lookup',
              algorithm: :concurrently,
              if_not_exists: true
  end
end
