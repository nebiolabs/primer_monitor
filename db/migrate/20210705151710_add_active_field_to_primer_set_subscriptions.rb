class AddActiveFieldToPrimerSetSubscriptions < ActiveRecord::Migration[6.1]
  def change
    add_column :primer_set_subscriptions, :active, :boolean, default: true, null: false
  end
end
