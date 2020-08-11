class CreateOligos < ActiveRecord::Migration[6.0]
  def change
    create_table :oligos do |t|
      t.string :name
      t.string :sequence
      t.references :amplicon, foreign_key: true
      t.ref_start :integer
      t.ref_end :integer
      
      t.timestamps
    end
  end
end
