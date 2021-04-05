class RemoveDuplicateLocationsAgain < ActiveRecord::Migration[6.1]
  def up
    execute "
    BEGIN TRANSACTION;

    with originals as (
    select count(*), min(id) as first_id, region, subregion, division, subdivision from detailed_geo_locations dgl 
    group by region, subregion, division, subdivision
    having count(id) > 1
    ),
    duplicates as (
    select id as duplicate_id, region, subregion, division, subdivision from detailed_geo_locations dgl 
    )
    update fasta_records
    set detailed_geo_location_id = cte.first_id
    from (
    select first_id, duplicates.duplicate_id from originals
    inner join duplicates on ((originals.region = duplicates.region) or (originals.region is null and duplicates.region is null))
               and ((originals.subregion = duplicates.subregion) or (originals.subregion is null and duplicates.subregion is null)) 
               and ((originals.division = duplicates.division) or (originals.division is null and duplicates.division is null))
               and ((originals.subdivision = duplicates.subdivision) or (originals.subdivision is null and duplicates.subdivision is null))
               and originals.first_id <> duplicates.duplicate_id
    inner join fasta_records on fasta_records.detailed_geo_location_id = duplicates.duplicate_id
    ) as cte
    where fasta_records.detailed_geo_location_id = cte.duplicate_id;
    
    with originals as (
    select count(*), min(id) as first_id, region, subregion, division, subdivision from detailed_geo_location_aliases dgl 
    group by region, subregion, division, subdivision
    having count(id) > 1
    ),
    duplicates as (
    select id as duplicate_id, region, subregion, division, subdivision from detailed_geo_location_aliases dgl 
    )
    update detailed_geo_locations
    set detailed_geo_location_alias_id = cte.first_id
    from (
    select first_id, duplicates.duplicate_id from originals
    inner join duplicates on ((originals.region = duplicates.region) or (originals.region is null and duplicates.region is null))
               and ((originals.subregion = duplicates.subregion) or (originals.subregion is null and duplicates.subregion is null)) 
               and ((originals.division = duplicates.division) or (originals.division is null and duplicates.division is null))
               and ((originals.subdivision = duplicates.subdivision) or (originals.subdivision is null and duplicates.subdivision is null))
               and originals.first_id <> duplicates.duplicate_id
    inner join detailed_geo_locations on detailed_geo_locations.detailed_geo_location_alias_id = duplicates.duplicate_id
    ) as cte
    where detailed_geo_locations.detailed_geo_location_alias_id = cte.duplicate_id;

    with originals as (
    select count(*), min(id) as first_id, region, subregion, division, subdivision from detailed_geo_locations dgl 
    group by region, subregion, division, subdivision
    having count(id) > 1
    ),
    duplicates as (
    select id as duplicate_id, region, subregion, division, subdivision from detailed_geo_locations dgl 
    )
    delete from detailed_geo_locations where id in (
    select duplicates.duplicate_id from originals
    inner join duplicates on ((originals.region = duplicates.region) or (originals.region is null and duplicates.region is null))
               and ((originals.subregion = duplicates.subregion) or (originals.subregion is null and duplicates.subregion is null)) 
               and ((originals.division = duplicates.division) or (originals.division is null and duplicates.division is null))
               and ((originals.subdivision = duplicates.subdivision) or (originals.subdivision is null and duplicates.subdivision is null))
               and originals.first_id <> duplicates.duplicate_id
     );
     
     with originals as (
      select count(*), min(id) as first_id, region, subregion, division, subdivision from detailed_geo_location_aliases dgl 
      group by region, subregion, division, subdivision
      having count(id) > 1
      ),
      duplicates as (
      select id as duplicate_id, region, subregion, division, subdivision from detailed_geo_location_aliases dgl 
      )
      delete from detailed_geo_location_aliases where id in (
      select duplicates.duplicate_id from originals
      inner join duplicates on ((originals.region = duplicates.region) or (originals.region is null and duplicates.region is null))
                 and ((originals.subregion = duplicates.subregion) or (originals.subregion is null and duplicates.subregion is null)) 
                 and ((originals.division = duplicates.division) or (originals.division is null and duplicates.division is null))
                 and ((originals.subdivision = duplicates.subdivision) or (originals.subdivision is null and duplicates.subdivision is null))
                 and originals.first_id <> duplicates.duplicate_id
       );

    COMMIT TRANSACTION;
    "
  end

  def down
    puts "Irreversible :("
  end
end
