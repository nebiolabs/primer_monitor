# frozen_string_literal: true

# variant_overlaps now keeps only recent sequences; this lets it find them, and max(date_collected), without
# scanning fasta_records. Built concurrently so the nightly loads aren't blocked on production.
class AddDateCollectedIndexToFastaRecords < ActiveRecord::Migration[8.1]
  disable_ddl_transaction!

  def change
    add_index :fasta_records, :date_collected, algorithm: :concurrently, if_not_exists: true
  end
end
