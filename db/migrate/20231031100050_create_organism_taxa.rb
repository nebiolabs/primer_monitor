class CreateOrganismTaxa < ActiveRecord::Migration[6.1]
  def change
    create_table :organism_taxa do |t|
      t.string :name, null: false, unique: true
      t.string :reference_accession, null: false
      t.references :organism, foreign_key: true, null:false
      t.string :caller_id, null: true
      t.timestamps
    end
  end
end
