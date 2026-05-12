class AddRefToVariantSites < ActiveRecord::Migration[7.1]
  def change
    add_column :variant_sites, :ref, :string
  end
end
