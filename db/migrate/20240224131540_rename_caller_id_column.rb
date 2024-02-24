class RenameCallerIdColumn < ActiveRecord::Migration[6.1]
  def change
    rename_column :organism_taxa, :caller_id, :lineage_caller_id
  end
end
