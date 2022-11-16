class AddSubmittedDateToFastaRecords < ActiveRecord::Migration[6.1]
  def self.up
    add_column :fasta_records, :submitted_date, :date
  end

  def self.down
    remove_column :fasta_records, :submitted_date
  end
end
