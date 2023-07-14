class RecreateLineages < ActiveRecord::Migration[6.1]
  def change

    drop_table :lineages

    create_table :lineages do |t|
      t.string :name, null: false, unique: true
      t.string :caller_name
      t.references :organism, foreign_key: true, null:false
      t.index :name, unique: true
      t.timestamps
    end
  end
end
