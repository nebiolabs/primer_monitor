class AddSubmittedDateToFastaRecords < ActiveRecord::Migration[6.1]
  def self.up
    add_column :fasta_records, :pangolin_version, :text
  end

  def self.down
    remove_column :fasta_records, :pangolin_version
  end
end
