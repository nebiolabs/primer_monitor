class GrantPermissionsOnNewTables < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      GRANT SELECT ON TABLE lineage_callers TO primer_monitor_ro;
      GRANT SELECT ON TABLE organism_taxa TO primer_monitor_ro;
      GRANT SELECT ON TABLE lineage_calls TO primer_monitor_ro;
    SQL
  end

  def down
    execute <<-SQL
      REVOKE SELECT ON TABLE lineage_callers FROM primer_monitor_ro;
      REVOKE SELECT ON TABLE organism_taxa FROM primer_monitor_ro;
      REVOKE SELECT ON TABLE lineage_calls FROM primer_monitor_ro;
    SQL
  end
end
