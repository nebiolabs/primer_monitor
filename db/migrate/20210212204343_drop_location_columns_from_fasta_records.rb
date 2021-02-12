class DropLocationColumnsFromFastaRecords < ActiveRecord::Migration[6.1]
  def change
    remove_column :fasta_records, :region, :string
    remove_column :fasta_records, :country, :string
    remove_column :fasta_records, :division, :string
  end
end
