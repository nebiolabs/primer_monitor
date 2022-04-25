class AddCreatedAtIndex < ActiveRecord::Migration[6.1]
  def change
    add_index :fasta_records, :created_at
  end
end
