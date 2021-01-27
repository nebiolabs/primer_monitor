class RenameAmpliconToPrimerSet < ActiveRecord::Migration[6.1]
  def change
    rename_table :amplicons, :primer_sets
    rename_column :oligos, :amplicon_id, :primer_set_id
  end
end
