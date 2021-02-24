class CreateGenomicFeatures < ActiveRecord::Migration[6.1]
  def change
    create_table :genomic_features do |t|
      t.string :name
      t.string :type
      t.integer :ref_start
      t.integer :ref_end
      t.references :organism, foreign_key: true, null:false, index: true

      t.timestamps
    end
  end
end
