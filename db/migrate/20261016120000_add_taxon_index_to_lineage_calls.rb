# frozen_string_literal: true

# Production already has index_lineage_calls_on_taxon, created by hand outside the migrations. Linking imported
# calls to fasta_records joins on taxon, so manage it here. Built concurrently so lineage loads aren't blocked.
class AddTaxonIndexToLineageCalls < ActiveRecord::Migration[8.1]
  disable_ddl_transaction!

  def change
    add_index :lineage_calls, :taxon, algorithm: :concurrently, if_not_exists: true
  end
end
