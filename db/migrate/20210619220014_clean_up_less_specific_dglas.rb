class CleanUpLessSpecificDglas < ActiveRecord::Migration[6.1]
  def up
    execute <<~SQL
    --correct fix for ireland #N 
    update fasta_records  set detailed_geo_location_id  = 4114 where detailed_geo_location_id = 16007;
    delete from detailed_geo_locations dgla2 where id = 16007;
    delete from detailed_geo_location_aliases dgla where region = '#N';
    
    --created too duplicate regions for these...  
    delete from detailed_geo_location_aliases dgla1
    where exists (select 1 from detailed_geo_location_aliases dgla2 where dgla2.id < dgla1.id and dgla1.region = dgla2.region
    and dgla1.subregion is null and dgla2.subregion is null);  
   
    delete from detailed_geo_location_aliases dgla1
    where exists (select 1 from detailed_geo_location_aliases dgla2 where dgla2.id < dgla1.id and dgla1.region = dgla2.region
    and dgla1.subregion = dgla2.subregion 
    and dgla1.division is null and dgla2.division is null );
   
    delete from detailed_geo_location_aliases dgla1
    where exists (select 1 from detailed_geo_location_aliases dgla2 
    where dgla2.id < dgla1.id and dgla1.region = dgla2.region
    and dgla1.subregion = dgla2.subregion 
    and dgla1.division = dgla2.division
    and dgla1.subdivision is null and dgla2.subdivision is null );

    SQL
  end
end
