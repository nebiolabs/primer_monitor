class AddShortNameToOligos < ActiveRecord::Migration[6.1]
  def change
    add_column :oligos, :short_name, :string
  end
end
