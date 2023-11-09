class AddTaxonIdToOligos < ActiveRecord::Migration[6.1]
  def change
    add_reference :oligos,  :organism_taxon, foreign_key: {to_table: :organism_taxa}, null:true
  end
end
