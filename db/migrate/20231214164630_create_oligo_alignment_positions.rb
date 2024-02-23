class CreateOligoAlignmentPositions < ActiveRecord::Migration[6.1]
  def up
    create_table :oligo_alignment_positions do |t|
      t.references :organism_taxon, foreign_key: {to_table: :organism_taxa}, null:false
      t.references :oligo, foreign_key: true, null:false
      t.integer :ref_start, null: false
      t.integer :ref_end, null: false
      t.timestamps
    end

    execute <<-SQL
        INSERT INTO oligo_alignment_positions (organism_taxon_id, oligo_id, ref_start, ref_end, created_at, updated_at)
        SELECT organism_taxon_id, id, ref_start, ref_end, NOW(), NOW() FROM oligos WHERE oligos.organism_taxon_id IS NOT NULL;
    SQL
  end

  def down
    raise ActiveRecord::IrreversibleMigration
  end
end
