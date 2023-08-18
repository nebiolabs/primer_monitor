class AddUniqueConstraintToPrimerSetName < ActiveRecord::Migration[6.1]
  def up
      add_index :primer_sets, [:name], unique: true
  end

  def down
    remove_index :primer_sets, column: [:name]
  end
end
