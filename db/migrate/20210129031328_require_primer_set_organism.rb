class RequirePrimerSetOrganism < ActiveRecord::Migration[6.1]
  def change
    change_column :primer_sets, :organism_id, :integer, :null => false
  end
end
