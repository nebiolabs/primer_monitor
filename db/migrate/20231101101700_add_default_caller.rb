class AddDefaultCaller < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      INSERT INTO lineage_callers(name, version_specifiers, script_name, created_at, updated_at) 
      VALUES ('default', NULL, 'default', NOW(), NOW());
    SQL
  end

  def down
    execute <<-SQL
      DELETE FROM lineage_callers WHERE name='default';
    SQL
  end
end
