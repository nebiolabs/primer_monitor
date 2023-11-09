class RenameOrganismTaxaIdToOrganismTaxonId < ActiveRecord::Migration[6.1]
  def change
    rename_column :variant_sites, :organism_taxa_id, :organism_taxon_id
    rename_column :fasta_records, :organism_taxa_id, :organism_taxon_id
  end
end