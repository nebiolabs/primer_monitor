class RenameSubmittedDateToDateSubmitted < ActiveRecord::Migration[6.1]
  def self.up
    rename_column :fasta_records, :submitted_date, :date_submitted
  end

  def self.down
    rename_column :fasta_records, :date_submitted, :submitted_date
  end
end