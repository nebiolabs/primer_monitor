class ChangeDateSubmittedName < ActiveRecord::Migration[6.1]
  def change
    rename_column :fasta_records, :date_submitted, :date_collected
  end
end
