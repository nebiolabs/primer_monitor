class AddLookbackToUsers < ActiveRecord::Migration[6.1]
    def change
      add_column :users, :lookback_days, :integer, null:false, default: 30
    end
end
