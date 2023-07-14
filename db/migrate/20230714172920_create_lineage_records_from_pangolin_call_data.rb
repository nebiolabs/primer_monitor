class CreateLineageRecordsFromPangolinCallData < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      INSERT INTO lineages (name, caller_name, organism_id, created_at, updated_at) 
      SELECT DISTINCT _lineage_name, 'pangolin', 1, NOW(), NOW() FROM pangolin_calls;

      UPDATE pangolin_calls SET lineage_id=(SELECT id FROM lineages WHERE name=pangolin_calls._lineage_name);
    SQL
    # (name, "pangolin", 1, NOW(), NOW())
  end

  def down
    raise ActiveRecord::IrreversibleMigration
  end

end
