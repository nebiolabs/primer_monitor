class AddUserIdColumn < ActiveRecord::Migration[6.1]
  def change
    add_reference :proposed_notifications, :user, null: false, foreign_key: true
  end
end
