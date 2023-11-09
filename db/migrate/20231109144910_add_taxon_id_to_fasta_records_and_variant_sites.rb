class AddTaxonIdToFastaRecordsAndVariantSites < ActiveRecord::Migration[6.1]
  def change
    remove_reference :fasta_records,  :organism, foreign_key: true, null:true
    remove_reference :variant_sites,  :organism, foreign_key: true, null:true

    add_reference :fasta_records,  :organism_taxa, foreign_key: true, null:true
    add_reference :variant_sites,  :organism_taxa, foreign_key: true, null:true
  end
end
d