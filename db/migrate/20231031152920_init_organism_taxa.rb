class InitOrganismTaxa < ActiveRecord::Migration[6.1]
  def up

    execute <<-SQL
      INSERT INTO organism_taxa(organism_id, name, reference_accession, ncbi_taxon_id, created_at, updated_at) 
      SELECT id, name, reference_accession, ncbi_taxon_id, NOW(), NOW() from organisms;
    SQL

    # reference accession belongs to organism_taxon now
    remove_column :organisms, :reference_accession
    remove_column :organisms, :ncbi_taxon_id

  end

  def down
    raise ActiveRecord::IrreversibleMigration
  end

end
