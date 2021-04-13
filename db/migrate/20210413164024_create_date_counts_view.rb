class CreateDateCountsView < ActiveRecord::Migration[6.1]
  def up
    execute "
    CREATE OR REPLACE VIEW date_counts AS
      SELECT COUNT(*) as total_count, detailed_geo_location_id, date_collected FROM fasta_records fr
      GROUP BY detailed_geo_location_id, date_collected
    "
  end

  def down
    execute <<-SQL
      DROP VIEW IF EXISTS date_counts;
    SQL
  end
end
