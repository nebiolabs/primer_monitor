class FixCallerIdColumn < ActiveRecord::Migration[6.1]
  def change
    remove_column :organism_taxa, :caller_id, :string
    add_reference :organism_taxa, :caller, foreign_key:  {to_table: :lineage_callers}, null:true
  end
end
