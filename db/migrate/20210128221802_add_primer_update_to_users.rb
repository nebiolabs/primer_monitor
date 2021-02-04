class AddPrimerUpdateToUsers < ActiveRecord::Migration[6.1]
  def change
    add_column :users, :send_primer_updates, :boolean, default: false, null: false
  end
end
