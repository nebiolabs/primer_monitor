class CreateProposedNotificationsTable < ActiveRecord::Migration[6.1]
  def change
    create_table :proposed_notifications do |t|
      t.references :primer_sets, foreign_key: true
      t.references :oligos, foreign_key: true
      t.references :verified_notifications, foreign_key: true
      t.integer :coordinate
      t.float :fraction_variant

      t.timestamps
    end
  end
end
