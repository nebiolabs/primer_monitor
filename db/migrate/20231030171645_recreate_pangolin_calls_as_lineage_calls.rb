class RecreatePangolinCallsAsLineageCalls < ActiveRecord::Migration[6.1]
  def change

    execute <<-SQL
      DROP MATERIALIZED VIEW lineage_info;
    SQL

    remove_reference :fasta_records, :pangolin_call, null: true, foreign_key: {to_table: :pangolin_calls}
    remove_reference :fasta_records, :pending_pangolin_call, null: true, foreign_key: {to_table: :pangolin_calls}

    drop_table :pangolin_calls

    create_table :lineage_calls do |t|
      t.string :taxon, null: false, unique: true
      t.references :lineages, foreign_key: true, null:true
      t.references :lineage_callers, foreign_key: true, null:true
      t.string :metadata, null: true
      t.timestamps
    end

    add_reference :fasta_records, :lineage_call, null: true, foreign_key: {to_table: :lineage_calls}
    add_reference :fasta_records, :pending_lineage_call, null: true, foreign_key: {to_table: :lineage_calls}
  end
end
