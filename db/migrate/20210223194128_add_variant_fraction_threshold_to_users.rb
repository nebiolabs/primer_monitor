class AddVariantFractionThresholdToUsers < ActiveRecord::Migration[6.1]
  def change
    add_column :users, :variant_fraction_threshold, :float, null: false, default: 0.1
  end
end
