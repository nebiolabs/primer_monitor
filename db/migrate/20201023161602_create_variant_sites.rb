class CreateVariantSites < ActiveRecord::Migration[6.0]
  def change
    create_table :variant_sites do |t|
      t.integer :position
      t.string :variant_type
      t.string :variant
      t.references :fasta_record, foreign_key: true

      t.timestamps
    end
  end
end
