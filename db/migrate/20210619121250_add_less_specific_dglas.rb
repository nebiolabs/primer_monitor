class AddLessSpecificDglas < ActiveRecord::Migration[6.1]
  def up
    execute <<~SQL
      -- convention is "USA"
      update detailed_geo_locations dgl 
      set detailed_geo_location_alias_id  = (select id from detailed_geo_location_aliases where division = dgla.division and dgla.subregion = 'USA')
      from detailed_geo_location_aliases dgla 
      where dgla.subregion = 'United States';
      delete from detailed_geo_location_aliases where subregion = 'United States';
      -- misspelled Eurpope
      update detailed_geo_locations dgl set detailed_geo_location_alias_id  = 5840 where dgl.detailed_geo_location_alias_id = 15748;
      delete from detailed_geo_location_aliases dgla2 where id = 15748;
      update detailed_geo_locations dgl set detailed_geo_location_alias_id  = 14067 where dgl.detailed_geo_location_alias_id = 15743;
      delete from detailed_geo_location_aliases dgla2 where id = 15743;
      update detailed_geo_locations dgl set detailed_geo_location_alias_id  = 14068 where dgl.detailed_geo_location_alias_id = 15747;
      delete from detailed_geo_location_aliases dgla2 where id = 15747;
      update detailed_geo_locations dgl set detailed_geo_location_alias_id  = 14070 where dgl.detailed_geo_location_alias_id = 15745;
      delete from detailed_geo_location_aliases dgla2 where id = 15745;
      update detailed_geo_locations dgl set detailed_geo_location_alias_id  = 17656 where dgl.detailed_geo_location_alias_id = 15749;
      delete from detailed_geo_location_aliases dgla2 where id = 15749;
      update detailed_geo_locations dgl set detailed_geo_location_alias_id  = 4732 where dgl.detailed_geo_location_alias_id = 15746;
      delete from detailed_geo_location_aliases dgla2 where id = 15746;
      update detailed_geo_locations dgl set detailed_geo_location_alias_id  = 5670 where dgl.detailed_geo_location_alias_id = 15744;
      delete from detailed_geo_location_aliases dgla2 where id = 15744;
      --#N record for ireland
      update detailed_geo_locations dgl set detailed_geo_location_alias_id  = 4114 where dgl.detailed_geo_location_alias_id = 16007;
      delete from detailed_geo_location_aliases dgla2 where id = 16007;
      -- South Americs
      update detailed_geo_locations dgl 
      set detailed_geo_location_alias_id  = (select id from detailed_geo_location_aliases 
      where region = 'South America' and subregion like '%Maarten' and division like '%Maarten' )
      where dgl.region = 'South Americs';
      delete from detailed_geo_location_aliases where region = 'South Americs';
      -- Euorpe
      update detailed_geo_locations dgl 
      set detailed_geo_location_alias_id  = (select id from detailed_geo_location_aliases 
      where region = 'Europe' and subregion = dgl.subregion and division = dgl.division and subdivision is null )
      where dgl.region = 'Euorpe';
      delete from detailed_geo_location_aliases where region = 'Euorpe';
      --header line
      delete from fasta_records where strain = 'strain';
      delete from detailed_geo_locations where region = 'region';
      delete from detailed_geo_location_aliases where region = 'region';
      -- Aisa
      update detailed_geo_locations dgl 
      set detailed_geo_location_alias_id  = (select id from detailed_geo_location_aliases 
      where region = 'Asia' and subregion = dgl.subregion and division = dgl.division and subdivision is null )
      where dgl.region = 'Aisa';
      delete from detailed_geo_location_aliases where region = 'Aisa';
      -- Switzerland not in Europe
      update detailed_geo_locations dgl 
      set detailed_geo_location_alias_id  = (select id from detailed_geo_location_aliases 
      where region = 'Europe' and subregion = dgl.region and division = dgl.division and subdivision is null )
      where dgl.region = 'Switzerland';
      delete from detailed_geo_location_aliases where region = 'Switzerland';

      --add less specific locations for easier subscription to an area of interest
      insert into detailed_geo_location_aliases (world, region, created_at, updated_at)
      select distinct world, region, now(), now() 
      from detailed_geo_location_aliases dgla
      where not exists (select 1 from detailed_geo_location_aliases dgla2 where dgla2.region=dgla.region and dgla.subregion is null);
      
      insert into detailed_geo_location_aliases (world, region, subregion, created_at, updated_at)
      select distinct world, region, subregion, now(),now() from detailed_geo_location_aliases dgla
      where not exists (select 1 from detailed_geo_location_aliases dgla2 where dgla2.region=dgla.region 
                and dgla2.subregion = dgla.subregion and dgla.division is null);
      
      insert into detailed_geo_location_aliases (world, region, subregion, division, created_at, updated_at)
      select distinct world, region, subregion, division, now(), now() from detailed_geo_location_aliases dgla
      where not exists (select 1 from detailed_geo_location_aliases dgla2 where dgla2.region=dgla.region 
                and dgla2.subregion = dgla.subregion and dgla2.division = dgla.division and dgla.subdivision is null);

    SQL
  end
end
