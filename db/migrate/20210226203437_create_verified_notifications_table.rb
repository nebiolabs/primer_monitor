class CreateVerifiedNotificationsTable < ActiveRecord::Migration[6.1]
  def change
    create_table :verified_notifications do |t|
      t.references :users, foreign_key: true
      t.date :date_sent
      t.string :status

      t.timestamps
    end
  end
end


