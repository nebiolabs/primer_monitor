class TransferLocationToGeoLocationTable < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      INSERT INTO geo_location(region, division, created_at, updated_at)
      SELECT DISTINCT region, SPLIT_PART(division, '-', 1), CURRENT_DATE, CURRENT_DATE FROM fasta_records;
      
      UPDATE fasta_records
      SET geo_loc_id = (SELECT id FROM geo_location WHERE region = fasta_records.region AND division = SPLIT_PART(fasta_records.division, '-', 1));
    
    SQL
  end
  def down
    puts "Irreversible :("
  end
end
