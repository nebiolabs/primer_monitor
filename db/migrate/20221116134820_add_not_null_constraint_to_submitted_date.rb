class AddNotNullConstraintToSubmittedDate < ActiveRecord::Migration[6.1]
  def self.up
    change_column_null :fasta_records, :submitted_date, false
  end

  def self.down
    change_column_null :fasta_records, :submitted_date, true
  end
end