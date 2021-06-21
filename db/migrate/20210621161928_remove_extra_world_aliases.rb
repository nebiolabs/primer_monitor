class RemoveExtraWorldAliases < ActiveRecord::Migration[6.1]
  def up
    execute <<~SQL
    --only need one "world"
    delete from detailed_geo_location_aliases dgla1
    where exists (select 1 from detailed_geo_location_aliases dgla2 where dgla2.id < dgla1.id and dgla1.world = dgla2.world
    and dgla1.region is null and dgla2.region is null);
    SQL
  end
end
