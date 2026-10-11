# frozen_string_literal: true

# users.active was never maintained (only seeds set it), so nearly every real user was false. It now marks
# whether the user gets email: inactive users can still sign in, but are sent no notifications.
# Every existing user starts active; deactivate departed users afterwards.
class MakeUsersActiveMeaningful < ActiveRecord::Migration[8.1]
  def up
    execute 'UPDATE users SET active = true'
    change_table :users, bulk: true do |t|
      t.change_default :active, from: false, to: true
      t.change_null :active, false
    end
  end

  def down
    change_table :users, bulk: true do |t|
      t.change_null :active, true
      t.change_default :active, from: true, to: false
    end
  end
end
