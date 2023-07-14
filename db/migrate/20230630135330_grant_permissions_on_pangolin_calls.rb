class GrantPermissionsOnPangolinCalls < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      GRANT SELECT ON TABLE pangolin_calls TO primer_monitor_ro;
    SQL
  end

  def down
    execute <<-SQL
      REVOKE SELECT ON TABLE pangolin_calls FROM primer_monitor_ro;
    SQL
  end
end
