class CreatePangolinCalls < ActiveRecord::Migration[6.1]
  def change
    create_table :pangolin_calls do |t|
      t.string :taxon, null: false
      t.string :lineage, null: false
      t.string :conflict
      t.float :ambiguity_score
      t.string :scorpio_call
      t.numeric :scorpio_support
      t.numeric :scorpio_conflict
      t.text :scorpio_notes
      t.string :version, null: false
      t.string :pangolin_version, null: false
      t.string :scorpio_version, null: false
      t.string :constellation_version, null: false
      t.string :is_designated, null: false
      t.string :qc_status
      t.string :qc_notes
      t.string :note

      t.timestamps
    end
  end
end
