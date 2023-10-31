class CreateLineageCallers < ActiveRecord::Migration[6.1]
  def change

    create_table :lineage_callers do |t|
      t.string :name, null: false, unique: true
      t.string :version_specifiers, null: true
      t.string :script_name, null: false
      t.timestamps
    end
  end
end
