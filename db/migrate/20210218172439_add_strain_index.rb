class AddStrainIndex < ActiveRecord::Migration[6.1]
  def change
    add_index :fasta_records, :strain, unique: true
  end
end
