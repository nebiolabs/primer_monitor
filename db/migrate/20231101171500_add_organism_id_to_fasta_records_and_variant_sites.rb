class AddOrganismIdToFastaRecordsAndVariantSites < ActiveRecord::Migration[6.1]
  def change
    add_reference :fasta_records,  :organism, foreign_key: true, null:true
    add_reference :variant_sites,  :organism, foreign_key: true, null:true
  end
end
