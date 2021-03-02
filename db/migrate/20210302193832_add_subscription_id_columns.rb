class AddSubscriptionIdColumns < ActiveRecord::Migration[6.1]
  def change
    change_column_null :proposed_notifications, :primer_sets_id, false
    change_column_null :proposed_notifications, :oligos_id, false
    change_column_null :proposed_notifications, :coordinate, false
    change_column_null :proposed_notifications, :fraction_variant, false

    add_reference :proposed_notifications, :subscribed_geo_locations, null: false, foreign_key: true
    add_reference :proposed_notifications, :primer_set_subscriptions, null: false, foreign_key: true
  end
end
