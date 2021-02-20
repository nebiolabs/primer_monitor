class AddAliasToOrganism < ActiveRecord::Migration[6.1]
  def change
    add_column :organisms, :alias, :string
  end
end
