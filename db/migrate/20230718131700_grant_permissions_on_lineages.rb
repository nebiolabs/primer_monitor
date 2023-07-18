class GrantPermissionsOnLineages < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      GRANT SELECT ON TABLE lineages TO primer_monitor_ro;
    SQL
  end

  def down
    execute <<-SQL
      REVOKE SELECT ON TABLE lineages FROM primer_monitor_ro;
    SQL
  end
end
