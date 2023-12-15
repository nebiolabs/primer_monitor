class RemovePositionColsFromOligos < ActiveRecord::Migration[6.1]
  def change
    remove_column :oligos, :ref_start, :integer
    remove_column :oligos, :ref_end, :integer
    remove_reference :oligos, :organism_taxon
  end
end
