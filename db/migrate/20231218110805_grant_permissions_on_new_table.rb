class GrantPermissionsOnNewTable < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      GRANT SELECT ON TABLE oligo_alignment_positions TO primer_monitor_ro;
    SQL
  end

  def down
    execute <<-SQL
      REVOKE SELECT ON TABLE oligo_alignment_positions FROM primer_monitor_ro;
    SQL
  end
end
