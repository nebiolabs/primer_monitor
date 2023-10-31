class AddTaxonIdToOrganismTaxa < ActiveRecord::Migration[6.1]
  def change
    add_column :organism_taxa, :ncbi_taxon_id, :integer
  end
end
