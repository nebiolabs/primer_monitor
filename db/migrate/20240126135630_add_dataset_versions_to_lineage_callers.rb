class AddDatasetVersionsToLineageCallers < ActiveRecord::Migration[6.1]
  def change
    add_column :lineage_callers,  :dataset_versions, :string
    add_column :lineage_callers,  :pending_dataset_versions, :string
  end
end