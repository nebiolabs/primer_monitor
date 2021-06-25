class OpenPermissionsOnViews < ActiveRecord::Migration[6.1]
  def up
    execute <<~SQL
      grant select on all tables in schema public to primer_monitor_ro
    SQL
  end
end
