class AddVariantBedLookbackToOrganisms < ActiveRecord::Migration[6.1]
  def change
    add_column :organisms,  :variant_bed_lookback_days, :integer, null: false, default: 180
  end
end