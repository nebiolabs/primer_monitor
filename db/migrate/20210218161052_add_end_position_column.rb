class AddEndPositionColumn < ActiveRecord::Migration[6.1]
  def change
    rename_column :variant_sites, :position, :ref_start
    add_column :variant_sites, :ref_end, :integer
  end
end
