class RenameLineagesIdToLineageId < ActiveRecord::Migration[6.1]
  def change
    rename_column :lineage_calls, :lineages_id, :lineage_id
  end
end