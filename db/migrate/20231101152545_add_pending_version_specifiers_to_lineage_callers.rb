class AddPendingVersionSpecifiersToLineageCallers < ActiveRecord::Migration[6.1]
  def change
    add_column :lineage_callers, :pending_version_specifiers, :string
  end
end
