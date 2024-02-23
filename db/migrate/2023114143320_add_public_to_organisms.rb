class AddPublicToOrganisms < ActiveRecord::Migration[6.1]
  def change
    add_column :organisms, :public, :boolean
  end
end
