class InitTaxonIdForOligos < ActiveRecord::Migration[6.1]
  def change
    execute <<-SQL
      UPDATE oligos SET organism_taxon_id=(SELECT id FROM organism_taxa WHERE ncbi_taxon_id=2697049) WHERE ref_start IS NOT NULL;
    SQL
  end
end
