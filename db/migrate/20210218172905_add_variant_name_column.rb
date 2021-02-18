class AddVariantNameColumn < ActiveRecord::Migration[6.1]
  def change
    add_column :fasta_records, :variant_name, :string
  end
end
