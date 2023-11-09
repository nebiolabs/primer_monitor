class FixLineageCallersIdColumn < ActiveRecord::Migration[6.1]
  def change
    rename_column :lineage_calls, :lineage_callers_id, :lineage_caller_id
  end
end
