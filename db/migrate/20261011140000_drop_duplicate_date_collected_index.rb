# frozen_string_literal: true

# Staging and production already had fasta_records_date_collected_idx, created by hand outside the migrations, so
# index_fasta_records_on_date_collected (20261011130000) duplicated it there. Drop the unmanaged copy; every insert
# into fasta_records was maintaining both.
class DropDuplicateDateCollectedIndex < ActiveRecord::Migration[8.1]
  disable_ddl_transaction!

  def up
    remove_index :fasta_records, name: 'fasta_records_date_collected_idx', algorithm: :concurrently, if_exists: true
  end

  def down
    # the managed index_fasta_records_on_date_collected covers the same column, so there is nothing to put back
  end
end
