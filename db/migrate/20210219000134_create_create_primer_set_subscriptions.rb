class CreateCreatePrimerSetSubscriptions < ActiveRecord::Migration[6.1]
  def change
    create_table :primer_set_subscriptions do |t|
      t.references :primer_set, foreign_key: true
      t.references :user, foreign_key: true

      t.timestamps
    end
    add_index :primer_set_subscriptions, [:user_id, :primer_set_id], unique: true
  end
end
